# Sleep Medicine Data Modeling and Analysis – Master's Thesis

This project presents an algorithm that automatically detects the drop in blood oxygen saturation (SpO₂) following each apnea event, determines the onset of that drop, and computes the time lag between apnea onset and desaturation onset. The main contribution is an automated method for measuring this lag on a large set of recordings from the Sleep Heart Health Study (SHHS) database.

## Method
1. **Reference value:** the local maximum of the SpO₂ signal within a 20-second window after apnea onset. 
2. **Desaturation onset:** within a 60-second window after apnea onset, the algorithm searches for a candidate point where SpO₂ falls below the reference value by at least 0.1%. The candidate is accepted only if the total drop reaches at least 3%. The algorithm then backtracks to the last point at the reference value before the drop, which is taken as the desaturation onset.
3. **Time lag:** the difference between desaturation onset and apnea onset (apnea onset is read from the XML annotation files).

## Sensitivity Analysis
To assess the robustness of the method, a sensitivity analysis was performed for the length of the reference-value window. The window was varied from 20 to 40 seconds (tested values: 20, 25, 30, 40) and the resulting mean time lag and number of detected episodes were compared. The mean lag changed by less than 0,5 seconds across the tested windows, indicating that the results are not strongly dependent on this parameter. 

## Results
The mean time lag between the cessation of breathing and the onset of SpO₂ desaturation, measured on 17,918 apnea episodes from 775 subjects, is 19.49 ± 4.85 s (mean ± SD).

<img width="865" height="751" alt="Distribution of time lag" src="https://github.com/user-attachments/assets/a79d8025-65dd-43c6-93f2-ef87275096de" />

## Technologies
Python, Jupyter, Pandas, NumPy, SciPy, bpyEDFlib

## Data
SHHS data are not included in this repository. They can be requested through the National Sleep Research Resource (NSRR): https://sleepdata.org/datasets/shhs
