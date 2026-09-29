"""
Mjerenje vremenskog kašnjenja između početka apneje i početka pada SpO2
na podacima iz baze SHHS-1 (Sleep Heart Health Study).

Podaci nisu uključeni u repozitorij. EDF i XML datoteke preuzete s platforme
NSRR (https://sleepdata.org/datasets/shhs) treba smjestiti u mapu data/.
"""

import os
import warnings
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pyedflib
from tqdm import tqdm

warnings.filterwarnings('ignore')

# ---- Putanje ----
BASE_PATH = os.path.join('data', 'polysomnography')
EDF_PATH = os.path.join(BASE_PATH, 'edfs', 'shhs1')
XML_PATH = os.path.join(BASE_PATH, 'annotations-events-nsrr', 'shhs1')

# ---- Parametri ----
LAG_MIN = 5.0       # epizode s kašnjenjem <= 5 s se odbacuju
REF_WINDOW_S = 20   # prozor za referentnu vrijednost (s); mijenja se za analizu osjetljivosti


def load_shhs_signals(edf_path):
    """Učitava SpO2 signal iz EDF datoteke; vrijednosti izvan 50-100 % postavlja na NaN."""
    try:
        f = pyedflib.EdfReader(edf_path)
        labels = [f.getLabel(i) for i in range(f.signals_in_file)]
        if 'SaO2' not in labels:
            f.close()
            return None, None
        idx = labels.index('SaO2')
        spo2 = f.readSignal(idx)
        fs = f.getSampleFrequency(idx)
        f.close()
        spo2 = np.where((spo2 >= 50) & (spo2 <= 100), spo2, np.nan)
        return np.round(spo2), fs
    except Exception:
        return None, None


def load_shhs_apneas(xml_path):
    """Čita apneje (bez hipopneja) iz XML anotacija: početak, trajanje i tip."""
    try:
        root = ET.parse(xml_path).getroot()
        apneas = []
        for e in root.findall('.//ScoredEvent'):
            ec = e.find('EventConcept')
            if ec is None:
                continue
            t = ec.text
            if 'Apnea' in t and 'Hypopnea' not in t and 'Arousal' not in t:
                start = float(e.find('Start').text)
                dur = float(e.find('Duration').text)
                if 'Obstructive' in t:
                    cls = 'OA'
                elif 'Central' in t:
                    cls = 'CA'
                elif 'Mixed' in t:
                    cls = 'MA'
                else:
                    cls = 'Unknown'
                apneas.append({'start': start, 'duration': dur, 'type': cls})
        return apneas
    except Exception:
        return []

def detect_desaturation_for_event(spo2_signal, fs, t0,
                                  ref_window_s=REF_WINDOW_S,
                                  max_onset_search=60, max_event_search=120,
                                  drop_threshold=3, onset_threshold_offset=0.1,
                                  recovery_pct=2.0, return_tolerance=0.5):
    """Detektira desaturaciju nakon apneje koja počinje u t0.
    Vraća početak pada, nadir, kraj desaturacije, referentnu vrijednost
    i veličinu pada, ili None ako desaturacija nije pronađena.
    """
    t0_idx = int(t0 * fs)
    if t0_idx >= len(spo2_signal):
        return None
    start = t0_idx
    end = min(len(spo2_signal), t0_idx + int(max_onset_search * fs))
    if end - start < 2:
        return None
    segment = spo2_signal[start:end]

    # 1) Referentna vrijednost: najviši SpO2 u prozoru ref_window_s, do prvog pada >= 3 %
    bw = min(int(ref_window_s * fs), len(segment))
    bseg = segment[:bw]
    run_max = np.nan
    run_max_rel = 0
    bi_rel = None
    for i in range(len(bseg)):
        v = bseg[i]
        if np.isnan(v):
            continue
        if np.isnan(run_max) or v > run_max:
            run_max = v
            run_max_rel = i
        elif run_max - v >= drop_threshold:
            bi_rel = run_max_rel
            break
    if bi_rel is None:
        # U prozoru nema pada >= 3 %: referenca je maksimum cijelog prozora
        if np.all(np.isnan(bseg)):
            return None
        run_max = float(np.nanmax(bseg))
        bi_rel = int(np.nanargmax(bseg))
    true_baseline = float(run_max)

    onset_threshold = true_baseline - onset_threshold_offset
    return_level = true_baseline - return_tolerance
    sab = segment[bi_rel:]     # pad se traži od referentne točke nadalje
    if len(sab) < 2:
        return None

    # 2) Kandidat za početak pada: prva točka ispod praga, uz ukupni pad >= 3 %
    drop_start = None
    search_from = 0
    while search_from < len(sab):
        cand = None
        for i in range(search_from, len(sab)):
            if sab[i] < onset_threshold:
                cand = i
                break
        if cand is None:
            break
        drop_ok = False
        m = sab[cand]
        for j in range(cand, len(sab)):
            vv = sab[j]
            if np.isnan(vv):
                continue
            if vv < m:
                m = vv
            if true_baseline - m >= drop_threshold:
                drop_ok = True
                break
            if vv >= return_level:
                break
        if drop_ok:
            drop_start = cand
            break
        search_from = cand + 1
    if drop_start is None:
        return None

    # 3) Povratak na posljednju točku s referentnom vrijednošću prije pada
    onset_rel = 0
    for j in range(drop_start - 1, -1, -1):
        if sab[j] >= true_baseline:
            onset_rel = j
            break
    onset_time = (start + bi_rel + onset_rel) / fs

    # 4) Nadir i kraj desaturacije (oporavak za recovery_pct iznad nadira)
    drop_start_abs = start + bi_rel + drop_start
    hard_limit = min(len(spo2_signal), drop_start_abs + int(max_event_search * fs))
    nadir_idx = drop_start_abs
    nadir_value = float(spo2_signal[drop_start_abs])
    offset_idx = hard_limit
    for k in range(drop_start_abs, hard_limit):
        vv = spo2_signal[k]
        if np.isnan(vv):
            continue
        if vv < nadir_value:
            nadir_value = vv
            nadir_idx = k
        elif vv >= nadir_value + recovery_pct:
            offset_idx = k
            break

    if true_baseline - nadir_value < drop_threshold:
        return None
    return {'onset_time': onset_time, 'nadir_time': nadir_idx / fs,
            'nadir_value': nadir_value, 'offset_time': offset_idx / fs,
            'baseline': true_baseline, 'drop': true_baseline - nadir_value}


def process_patient(patient_id, ref_window_s=REF_WINDOW_S):
    """Obrađuje jednog pacijenta i vraća kašnjenja za sve pronađene epizode."""
    edf_file = os.path.join(EDF_PATH, f'shhs1-{patient_id}.edf')
    xml_file = os.path.join(XML_PATH, f'shhs1-{patient_id}-nsrr.xml')
    if not (os.path.exists(edf_file) and os.path.exists(xml_file)):
        return None
    spo2, fs = load_shhs_signals(edf_file)
    if spo2 is None:
        return None
    apneas = sorted(load_shhs_apneas(xml_file), key=lambda a: a['start'])
    if len(apneas) == 0:
        return {'patient_id': patient_id, 'n_apneas': 0, 'n_matched': 0, 'lag_times': []}

    lag_results = []
    for apnea in apneas:
        t0 = apnea['start']
        d = detect_desaturation_for_event(spo2, fs, t0, ref_window_s=ref_window_s)
        if d is None:
            continue
        lag = d['onset_time'] - t0
        if lag > LAG_MIN:
            lag_results.append({
                'patient_id': patient_id, 'apnea_start': t0,
                'apnea_duration': apnea['duration'], 'apnea_type': apnea['type'],
                'lag_time': lag, 'baseline': d['baseline'], 'nadir': d['nadir_value'],
                'nadir_time': d['nadir_time'], 'offset_time': d['offset_time'],
                'drop': d['drop']})
    return {'patient_id': patient_id, 'n_apneas': len(apneas),
            'n_matched': len(lag_results), 'lag_times': lag_results}


if __name__ == '__main__':
    print(f"SHHS – vremensko kašnjenje (referentni prozor = {REF_WINDOW_S} s)")

    patient_ids = [f'{200000 + i:06d}' for i in range(1, 1001)]
    all_results, all_lag_times = [], []
    for pid in tqdm(patient_ids, desc="Obrada"):
        r = process_patient(pid)
        if r is not None:
            all_results.append(r)
            all_lag_times.extend(r['lag_times'])

    # Rezultati po epizodi
    df_all = pd.DataFrame(all_lag_times)
    df_all.to_csv('shhs_lag_events.csv', index=False)
    print(f"Spremljeno epizoda: {len(df_all)}")

    # Sažetak po pacijentu
    summary = []
    for r in all_results:
        if r['n_matched'] > 0:
            lags = [lt['lag_time'] for lt in r['lag_times']]
            summary.append({'patient_id': r['patient_id'], 'n_apneas': r['n_apneas'],
                            'n_matched': r['n_matched'], 'mean_lag': np.mean(lags),
                            'std_lag': np.std(lags), 'median_lag': np.median(lags)})
    df_summary = pd.DataFrame(summary)
    df_summary.to_csv('shhs_lag_summary.csv', index=False)
    print(f"Spremljeno pacijenata: {len(df_summary)}")
    print(f"Po epizodama: mean {df_all['lag_time'].mean():.2f} s, "
          f"median {df_all['lag_time'].median():.2f} s")
    print(f"Po pacijentu: mean {df_summary['mean_lag'].mean():.2f} s, "
          f"median {df_summary['mean_lag'].median():.2f} s, n = {len(df_summary)}")