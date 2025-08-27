def simulate_radar_measurements(state, radars):
RADARS = [
{
"id": "dc_east_10km",
"lat_deg": 38.907189, # computed
"lon_deg": -76.923629, # computed
"alt_m": 0.0,
"update_rate_hz": 2.0,
"sigma_range_m": 50.0,
"sigma_az_rad": 0.0009, # ~0.05°
"sigma_el_rad": 0.0009,
"sigma_doppler_mps": 3.0,
"dropout_prob": 0.05,
},
{
"id": "nyc_west_10km",
"lat_deg": 40.730602, # computed
"lon_deg": -74.046727, # computed
"alt_m": 0.0,
"update_rate_hz": 2.0,
"sigma_range_m": 50.0,
"sigma_az_rad": 0.0009,
"sigma_el_rad": 0.0009,
"sigma_doppler_mps": 3.0,
"dropout_prob": 0.05,
},
{
"id": "midpoint_west_15km",
"lat_deg": 39.825905, # computed
"lon_deg": -75.694756, # computed
"alt_m": 0.0,
"update_rate_hz": 2.0,
"sigma_range_m": 50.0,
"sigma_az_rad": 0.0009,
"sigma_el_rad": 0.0009,
"sigma_doppler_mps": 3.0,
"dropout_prob": 0.05,
},
{
"id": "pennslyvania_radar",
"lat_deg": 40.051626, # supplied
"lon_deg": -77.344419, # supplied
"alt_m": 0.0,
"update_rate_hz": 2.0,
"sigma_range_m": 50.0,
"sigma_az_rad": 0.0009,
"sigma_el_rad": 0.0009,
"sigma_doppler_mps": 3.0,
"dropout_prob": 0.05,
},
]
