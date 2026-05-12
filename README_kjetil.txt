#Collection of scripts and recipies for preparing land use data for NorESM-FATES simulations

Tools and data can be found on nac/qbo:
Tools: /div/no-backup-nac/users/kjetisaa/
Data: /div/no-backup-nac/users/LUH2/

Relevant weight files:
ne16: /div/no-backup-nac/users/masan/NorESM_stuff/mapping_weights_etc/map_0.5x0.5_GRDC_to_ne16np4_nomask_aave_da_test251104_3.nc
ne30: /div/no-backup-nac/users/masan/NorESM_stuff/mapping_weights_etc/map_fv0.9x1.25_to_ne30pg3_scripgrids_nomask_conserve_c251203.nc

Intermediate regridding files:
ne16: /div/no-backup-nac/users/masan/NorESM_stuff/mapping_weights_etc/griddata_360x720_070122.nc
ne30: /div/no-backup-nac/users/masan/NorESM_stuff/mapping_weights_etc/griddata_0.9x1.25_c070928.nc

Concrete steps (change to correct resolution and input files etc):
Step 1:
- copy surface data from Betzy to qbo. 
- download LUH2 and CLM5 THESIS datasets (following documentation). 
- load xesmf environment: 'module load Anaconda3/2022.05', 'source activate /div/no-backup-nac/conda_env/xesmf_env/'
- run tools-fates-landusedata (from src/): 
    - python -m landusedata lupft ../CLM_surfdata/surfdata_4x5_hist_1850_16pfts_c241007.nc ../LUH2/staticData_quarterdeg.nc ../LUH2/CLM5_current_luhforest_deg025.nc ../LUH2/CLM5_current_luhpasture_deg025.nc ../LUH2/CLM5_current_luhother_deg025.nc ../LUH2/CLM5_current_surf_deg025.nc
    - python -m landusedata luh2 ../CLM_surfdata/surfdata_4x5_hist_1850_16pfts_c241007.nc ../LUH2/staticData_quarterdeg.nc ../LUH2/states.nc ../LUH2/transitions.nc ../LUH2/management.nc
- change format: nccopy -k cdf5 fates_landuse_pft_map_4x5.nc fates_landuse_pft_map_4x5_cdf5.nc
- run constantlanduse:
    - python calculate_luh2_secondary_agedist.py ../tools-fates-landusedata/src/"file".nc
- change format: nccopy -k cdf5 fates_landuse_pft_map_4x5.nc fates_landuse_pft_map_4x5_cdf5.nc
