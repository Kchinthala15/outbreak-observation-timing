# Discrete shifts in outbreak response intensity

Reproducible change-point analysis of onset-to-hospitalisation delay during the
2018–2020 Ebola outbreak in the Democratic Republic of the Congo.

## What this tests

Surveillance systems do not only observe epidemics; they change how epidemics are
observed. When clinical suspicion forms, reporting and response intensity can shift
sharply rather than gradually — and if that shift is discrete, early-warning statistics
computed over fixed windows of observations lose resolution precisely during the
pre-recognition period when detection matters most.

This repository tests one specific, falsifiable part of that argument at the level of
individual patients: **within an ongoing outbreak, does response speed improve in
discrete steps, or as a smooth ramp?**

The measure is onset-to-hospitalisation delay — the interval between a patient
developing symptoms and reaching hospital.

## Method

An exhaustive change-point scan over every observed onset time that leaves at least
20 hospitalised cases on each side of the split. A two-sided Mann–Whitney U test is
computed at each candidate, and the Bonferroni correction uses the number of candidates
actually evaluated (not a nominal figure).

Whether the resulting pattern is stepwise or smooth is then adjudicated by OLS model
comparison — a step model against a linear trend, by AIC — plus a residual trend test
to confirm nothing gradual remains once the steps are accounted for.

## Findings

Wave 3 (n = 748 hospitalised cases, 709 candidate split points):

| Break | Mean delay before | Mean delay after | Bonferroni-corrected p |
|---|---|---|---|
| Day 77.8 | 4.88 d (n = 610) | 3.42 d (n = 138) | ≈ 0.002 |
| Day 12.0 | 6.08 d (n = 85) | 4.42 d (n = 663) | ≈ 0.016 |

Both breaks run in the predicted direction. A ±10-day window around the earlier break
shows the same contrast (6.00 d vs 4.60 d), ruling out a small-early-sample artefact.

Stepwise or smooth:

| Model | R² | AIC |
|---|---|---|
| Linear trend only | 0.019 | 3957.3 |
| Two steps only | 0.043 | 3940.4 |
| Two steps + linear | 0.047 | 3939.9 |

The two-step model outperforms a linear trend by roughly 17 AIC. Adding a linear term
to the steps yields no meaningful improvement, and after both steps are fitted no
residual time trend remains (Spearman ρ = 0.03, p = 0.45).

Waves 1 and 2, analysed under identical methodology, show no significant break
(corrected p = 1.00 each) — evidence that the Wave 3 result is a wave-specific signal
rather than a tendency of the search to find breaks wherever it looks.

**Interpretation.** Response intensity in Wave 3 improved through two discrete
transitions, not a continuous ramp.

## Data

The wave-level data are not redistributed here. Download them from Zenodo:

> Zenodo DOI [10.5281/zenodo.6104188](https://doi.org/10.5281/zenodo.6104188)

You need `EVD_Wave1.csv`, `EVD_Wave2.csv` and `EVD_Wave3.csv`. Each has columns
`onset`, `hosp`, `event1`, `event2`, with times in days from the start of the wave.
`event2 == 1` marks a hospitalised case; delay is computed as `hosp - onset`.

Place them in a directory of your choice and pass that path with `--data-dir`.

## Running it

```bash
pip install -r requirements.txt
python level4_changepoint_analysis.py --data-dir ./data
```

Options:

- `--data-dir` — directory holding the three wave CSVs (default: current directory)
- `--min-group` — minimum observations either side of a candidate split (default: 20)

The `--min-group` setting changes how many candidate splits are admissible, which in
turn changes the Bonferroni denominator. It is exposed deliberately so the correction
is transparent rather than buried.

## Reproducibility note

Every figure quoted above is printed by the script. Results were produced with
pandas 3.0.2, numpy 2.4.4, scipy 1.17.1 and statsmodels 0.14.6 under Python 3.

## Source data citation

Vossler H, et al. *Scientific Reports*, 2022. Data: Zenodo
[10.5281/zenodo.6104188](https://doi.org/10.5281/zenodo.6104188).

Please cite the original authors when using these data. This repository contains
analysis code only.

## Licence

MIT — see [LICENSE](LICENSE).
