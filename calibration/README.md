# Camera calibrations

`calibration_history.csv` lists every extrinsic calibration (EO) of the two
Marconi cameras, from the station's calibration archive, with the camera
position and angles (CIRN convention: azimuth clockwise from grid north,
tilt 0 = nadir / 90 = horizon, roll positive counter-clockwise seen from
behind). Positions are UTM 19N, elevations NAVD88.

Setups (the station ID in the image file names):

| Setup | Calibrations | Notes |
|---|---|---|
| CACO02 | 2021-03, 2023-03, 2024-01, 2024-06 | ~17 m from the later positions; lens IO 2021-02-25, 2024-01-05 |
| CACO03 | 2024-08-27, 2024-10-23, **2025-01-23** | lens IO 2024-08-01 (the same IO is used for every later setup). Both cameras turned ~13 deg between the October and January calibrations. 2025-01-23 is the day of the Jan 2025 lidar; that calibration was solved from ground control surveyed that day, NOT from the lidar, so January frames placed with it can be tested against the lidar independently. |
| CACO04 | 2025-02-19, 2025-11-04 | set up 2025-01-24: cameras moved ~5.7 m and turned back ~15 deg. `CACO05_<cam>_20250219_EO.yaml` holds the same values under an older name. |
| CACO05 | 2025-11-13 (`_EO-CV`) | the live station EO. c1 has the 2025-11-04 values; c2 moved 3.1 m. |

`c1_geom_sigma.csv`, `c2_geom_sigma.csv`: how far off the live CACO05
2025-11-13 calibration can place each pixel on the ground, written by
`georectify_uncertainty.py` from that solve's own targets
(`2025-11-13_Marconi_Extrinsic_Targets_<cam>_xyz.csv`, `gcp_uv_<cam>.csv`):
the solve repeated with each target left out (in-sample RMS 0.42 m for
both cameras; left out, c1 0.59 m, c2 0.41 m), then its parameter
uncertainty drawn 2000 times through `georectify.pixel_to_ground`. One row
per pixel column (every 256 px), ground range (20-1000 m) and plane
(z -1, 0, 1, 2 m): sigma of the easting and northing, cross-shore sigma
and the cross-shore pixel footprint; the `#` header holds the provenance
(inputs by sha256, the per-target leave-one-out errors, the covariance).
Read only by `dem_from_contours.py --sigma-obs` (opt-in). For the
2025-11-13 calibration only: rerun the script after any new calibration.

For historical processing (process_chelsea.py), point the rows of
`chelsea_setups.csv` at the calibration of the period: e.g. Jan 18-23 2025
-> `CACO03_<cam>_20250123_EO.yaml`. A survey-fitted or carried pointing
(fit_eo_to_survey.py, apply_pointing_correction.py) keeps the camera
position fixed, so after a move it absorbs the move into the angles
(the Jan 2025 lidar fit wanted ~22 deg of pan where the calibration has
~15 deg and a 5.8 m move); use it only where no calibration exists.

Timex bias corrections of the detector, one file per water-level timing mode
(`timex_timing.py`, `TIMING_MODE`; docs/VALIDATION.md section 2.7):

| File | Mode | What |
|---|---|---|
| `timex_bias_legacy.json` | `legacy` (the default) | The detector's built-in v6.8 values (c1 constant -11.03 px, c2 a 10-point curve), fitted with the temporal blend on. Never edited: the legacy products are the v6.8 ones, and `test_timing_mode.py` fails on any change. |
| `timex_bias_mid.json` | `mid` | Not shipped. Written only by `recalibrate_timing.py --steps derive` on the NUC, from the hand-clicked ground truth detected with the blend off; a previous one is kept as `.bak`. Without it a mid timex run is refused. |
