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
| CACO03 | 2024-08-27, 2024-10-23, **2025-01-23** | lens IO 2024-08-01 (the same IO is used for every later setup). Both cameras turned ~13 deg between the October and January calibrations. 2025-01-23 is the day of the Jan 2025 lidar. |
| CACO04 | 2025-02-19, 2025-11-04 | set up 2025-01-24: cameras moved ~5.7 m and turned back ~15 deg. `CACO05_<cam>_20250219_EO.yaml` holds the same values under an older name. |
| CACO05 | 2025-11-13 (`_EO-CV`) | the live station EO. c1 has the 2025-11-04 values; c2 moved 3.1 m. |

For historical processing (process_chelsea.py), point the rows of
`chelsea_setups.csv` at the calibration of the period: e.g. Jan 18-23 2025
-> `CACO03_<cam>_20250123_EO.yaml`. A survey-fitted or carried pointing
(fit_eo_to_survey.py, apply_pointing_correction.py) keeps the camera
position fixed, so after a move it absorbs the move into the angles
(the Jan 2025 lidar fit wanted ~22 deg of pan where the calibration has
~15 deg and a 5.8 m move); use it only where no calibration exists.
