#######################################
### Cyclone Harold CLIMADA analysis ###
#######################################

#%matplotlib inline
import numpy as np
from climada.hazard import Centroids, TropCyclone, TCTracks
import numpy as np
from climada.entity import LitPop
from climada.engine import Impact
from climada.entity import ImpactFuncSet, ImpfTropCyclone

# Import the tracks
tr_harold = TCTracks.from_ibtracs_netcdf(provider='usa', storm_id='2020092S09155')

# Plot the category progression
ax = tr_harold.plot();
ax.set_title('Cyclone Harold Category Progression', fontsize = 20, pad = 20) 

# Wind Speed Visualisation
harold_track = tr_harold.get_track('2020092S09155') 
lines = harold_track.max_sustained_wind.plot();

# Generate probabilistic events
tr_harold.equal_timestep()
tr_harold.calc_perturbed_trajectories(nb_synth_tracks=5)
tr_harold.plot();

# Select a specific track
tc_syn = tr_harold.get_track('2020092S09155_gen1')

# Compute the maximum sustained wind for each day.
print('Daily max sustained wind:', tc_syn.max_sustained_wind.groupby('time.day').max())

# Construct the centroids
min_lat, max_lat, min_lon, max_lon = -17.7, -10.2, 159, 174.9
cent = Centroids.from_pnt_bounds((min_lon, min_lat, max_lon, max_lat), res=0.5)
cent.check()

# Construct Cyclone Harold
tc_harold = TropCyclone.from_tracks(tr_harold, centroids=cent)
tc_harold.check()

# Plot the intensity map
ax = tc_harold.plot_intensity('2020092S09155');  # Harold
ax.set_title('Cyclone Harold Intensity Map', fontsize = 20, pad = 20)

#
### Video of Cyclone Harold - uncomment to run
#

# WARNING: executing the below will fail unless there is enough memory available (> 10G)
 
track_name = '2020092S09155' 
tr_harold = TCTracks.from_ibtracs_netcdf(provider='usa', storm_id='2020092S09155')
 
lon_min, lat_min, lon_max, lat_max = 159, -17.7, 174.9, -10.2
centr_video = Centroids.from_pnt_bounds((lon_min, lat_min, lon_max, lat_max), 0.04)
centr_video.check()
 
tc_video = TropCyclone()
 
tc_list, tr_coord = tc_video.video_intensity(track_name, tr_harold, centr_video, file_name='/Users/james/climada/results/harold_tc_fl.gif')

#
### Impact calculations
#

# VUT, FJI, SLB, TON with resolution 1km and financial_mode = income group.
exp_lp = LitPop.from_countries(countries=['VUT', 'FJI', 'SLB', 'TON'], res_arcsec=30, fin_mode='income_group')
exp_lp.check()
exp_lp.gdf.head()

# Visualising the exposure on a map
exp_lp.plot_raster()
tr_harold.equal_timestep(0.5)

# Calculate the centroids according to the exposures position
lat = exp_lp.gdf['latitude'].values
lon = exp_lp.gdf['longitude'].values
centrs = Centroids.from_lat_lon(lat, lon)
centrs.check()

tc_harold = TropCyclone.from_tracks(tr_harold, centroids=centrs)
tc_harold.check()

# Impact function for tropical cyclones
impf_tc = ImpfTropCyclone.from_emanuel_usa()

# Adding the impact function to an Impact function set
impf_set = ImpactFuncSet()
impf_set.append(impf_tc)
impf_set.check()

# Extracting the hazard type and id
[haz_type] = impf_set.get_hazard_types()
[haz_id] = impf_set.get_ids()[haz_type]

# Rename column and allocate id
exp_lp.gdf.rename(columns={"impf_": "impf_" + haz_type}, inplace=True)
exp_lp.gdf['impf_' + haz_type] = haz_id
exp_lp.check()
exp_lp.gdf.head()

# Calculate impact
imp = Impact()
imp.calc(exp_lp, impf_set, tc_harold, save_mat=True) 

print(f"The aggregated average annual impact: {round(imp.aai_agg,0)} $")
imp.at_event
print(imp.imp_mat)
imp.tot_value
imp.imp_mat()
imp.plot_hexbin_eai_exposure(buffer=1);
imp.plot_raster_eai_exposure()
imp.plot_rp_imp(buffer = 200)
imp.plot_scatter_eai_exposure()

