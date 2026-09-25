# Optical wave gauge, CACO05 camera 2, Dec 2024 - Mar 2025

Leak-free evaluation: split by whole days, validation never resampled,
quality-filtered, ADCP in view of the camera.

| experiment            | RMSE (m) | R2   |
|-----------------------|----------|------|
| timex                 | 0.340    | 0.68 |
| timex + water level   | 0.368    | 0.64 |
| bright                | 0.296    | 0.70 |
| bright + ImageNet     | 0.299    | 0.69 |
| bright + glare filter | 0.291    | 0.71 |
| + Inception-ResNetV2  | 0.288    | 0.71 |

Data changes: 0.340 -> 0.291. Model changes: 0.291 -> 0.288.

## Temporal averaging
Error autocorrelation between consecutive frames +0.74, so most of the
remaining error is shared within a day rather than independent noise.
Comparing like with like (same rows, averaged vs single-frame):
3-frame 0.252 -> 0.232 (8%); 7-frame 0.201 -> 0.178 (11%).
A naive comparison against the 0.288 headline suggests 38%, but that
mostly reflects the averaging windows keeping only long unbroken
daytime runs, which are the easier frames.

## Limits
Image sharpness falls 20-fold from calm to storm (median 122 at
<0.5 m, 6 at >2.5 m): spray and haze obscure the camera exactly when
waves are largest. Wave period was not predictable (R2 0.10); peak
period changes 1.70 s RMS per hour against a total spread of 2.45 s,
so a single time-averaged image cannot resolve it.

## What quality filtering is worth
Running the trained model over frames it never saw, split by whether
they passed quality control:

| frames                | n   | RMSE (m) | bias (m) |
|-----------------------|-----|----------|----------|
| passed QC (validation)| 360 | 0.288    | -0.054   |
| rejected by QC        | 954 | 1.974    | +1.398   |

Rejected frames are dark, hazy or glare-blown. Normalised to zero mean
and unit variance, a dark frame becomes amplified noise, and the model
reads that texture as heavy breaking: it over-predicts by 1.4 m on
average and fails confidently rather than visibly. Without the filter
the gauge is not usable; with it, it measures.
