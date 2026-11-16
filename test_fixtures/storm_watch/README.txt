STORM WATCH FORECAST FIXTURES: ALL SYNTHETIC (NOT REAL RESPONSES)

Read by test_storm_watch.py (T2). Written in Nov 2026 from the providers' documentation and from
memory of the formats, because no request to the forecast hosts could be made from the
development sandbox (its network policy refused them). The values are invented: a north-east
gale building to ~5 m offshore on 18 Nov 2026. They test the parsers' handling of the documented
shapes; they are not data, and they do not show that the real responses look like this.

  nws_points_SYNTHETIC.json   GET https://api.weather.gov/points/41.89,-69.96
                              gridX/gridY are placeholders (101,101): storm_watch.py never
                              hardcodes them; it caches the forecastGridData URL weekly.
  nws_grid_SYNTHETIC.json     GET <forecastGridData>: layers {uom, values:[{validTime, value}]}.
                              The wavePeriod unit and whether the beach cell has wave values
                              are NOT verified.
  nws_alerts_SYNTHETIC.json   GET https://api.weather.gov/alerts/active?zone=ANZ254,MAZ022
                              a Gale Warning (ANZ254), a Coastal Flood Watch (MAZ022), a cancelled
                              High Surf Advisory and a Test message.
  om_marine_SYNTHETIC.json    GET marine-api.open-meteo.com/v1/marine?latitude=41.90,40.50&...
                              two points -> a JSON array. 3-hourly here (the real one is hourly).
  om_wind_SYNTHETIC.json      GET api.open-meteo.com/v1/forecast?...&wind_speed_unit=ms
  gfswave_44008_SYNTHETIC.bull  GET polar.ncep.noaa.gov/waves/WEB/gfswave.latest_run/plots/gfswave.44008.bull
                              layout reconstructed from memory; column widths not verified (the
                              parser splits on '|' and whitespace, never on fixed columns).

Before trusting the parsers, capture one real response of each ON THE NUC:
    ./storm_watch.sh --dry-run --save-raw /tmp/sw_raw
    python3 test_storm_watch.py --real-captures /tmp/sw_raw
