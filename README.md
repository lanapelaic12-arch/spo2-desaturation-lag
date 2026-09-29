# Modeliranje i analiza podataka u medicini spavanja – diplomski rad

Razvijen je algoritam koji iz signala SpO₂ (zasićenost krvi kisikom) automatski prepoznaje pad nakon svake apneje i određuje trenutak njegova početka, a zatim računa vremensko kašnjenje kao razliku između početka pada i početka apneje. Glavni doprinos rada je automatizirana metoda mjerenja kašnjenja na velikom skupu zapisa iz baze SHHS (Sleep Heart Health Study).

## Način rada
1. Vrijednosti izvan 50-100% se odbacuju i zamijenjuju Nan vrijednostima.

2. **Referentna vrijednost:** lokalni maksimum signala SpO₂ unutar prozora od 20 sekundi nakon početka apneje.
3. **Početak pada:** unutar prozora od 60 sekundi nakon početka apneje traži se kandidat, točka u kojoj SpO₂ padne ispod referentne vrijednosti za najmanje 0,1 %. Kandidat se prihvaća ako ukupni pad dosegne najmanje 3 %. Algoritam se tada vraća na posljednju točku s referentnom vrijednošću prije pada i ta točka postaje početak pada.
4. **Vremensko kašnjenje:** razlika između početka pada SpO₂ i početka apneje (početak apneje čita se iz XML datoteka s anotacijama).

## Analiza osjetljivosti
Kako bi se provjerila robusnost metode, provedena je analiza osjetljivosti na duljinu prozora za određivanje referentne vrijednosti. Testirani su prozori od 20, 25, 30 i 40 sekundi te su uspoređeni dobivena srednja vrijednost kašnjenja i broj detektiranih epizoda. Srednje kašnjenje promijenilo se za manje od 0,5 sekundi, što pokazuje da rezultati ne ovise o ovom parametru. 

## Rezultati
Srednja vrijednost vremenskog kašnjenja između prestanka disanja i početka pada SpO₂, izmjerena na 17 918 apnejskih epizoda kod 775 ispitanika, iznosi 19,49 ± 4,85 s (srednja vrijednost ± standardna devijacija).

<img width="865" height="751" alt="Raspodjela vremenskog kašnjenja" src="https://github.com/user-attachments/assets/a79d8025-65dd-43c6-93f2-ef87275096de" />

## Tehnologije
Python, Jupyter, Pandas, NumPy, SciPy, pyEDFlib

## Podaci
Podaci iz baze SHHS nisu uključeni u repozitorij, a mogu se zatražiti putem platforme NSRR: https://sleepdata.org/datasets/shhs
