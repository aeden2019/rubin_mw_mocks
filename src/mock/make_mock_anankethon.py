"""
Create ANANKE mock catalog from FIRE data with the intent of replicating the anankethon results
(see https://github.com/athob/anankethon/blob/main/usage_example.ipynb)

"""

import numpy as np
import args
import gizmo_analysis as ga
import allvariables
from astropy import units as u
import ananke as an
import plot_ananke_inputs



if __name__ == "__main__":
 
    # Get the parameters from the config file
    params = allvariables.getparams()

    # Define parameter names
    param_names = [
        "snap", "sim", "sim_dir", "sat_id_dir", "rmin", "rmax", 
        "bmin", "bmax", "overdensity", "sat_mask", "only_sat_mask", 
        "halo_mask", "rand_mask", "subsample", "ananke_name", 
        "photo_sys", "cmd_magnames", "app_mag_lim_lo", "app_mag_lim_hi", 
        "abs_mag_lim_lo", "abs_mag_lim_hi", "ananke_r_max"
    ]

    # Initialize parameter dictionary
    params_dict = {}

    # Assign values to parameters and declare them locally
    for param_name, param_value in zip(param_names, params):
        params_dict[param_name] = param_value
        locals()[param_name] = param_value
    
#     # Declare all the parameters locally
#     snap = params[0]
#     sim = params[1]
#     sim_dir = params[2]
#     sat_id_dir = params[3]
#     rmin = params[4]
#     rmax = params[5]
#     bmin = params[6]
#     bmax = params[7]
#     overdensity = params[8]
#     sat_mask = params[9]
#     only_sat_mask = params[10]
#     halo_mask = params[11]
#     rand_mask = params[12]
#     subsample = params[13]
#     ananke_name = params[14]
#     photo_sys = params[15]
#     cmd_magnames = params[16]
#     app_mag_lim_lo = params[17]
#     app_mag_lim_hi = params[18]
#     abs_mag_lim_lo = params[19]
#     abs_mag_lim_hi = params[20]
#     ananke_r_max = params[21]
    
    # Print parameters (for debugging)
    print(f"Parameters:")
    print(f"  snap: {snap}")
    print(f"  sim: {sim}")
    print(f"  sim_dir: {sim_dir}")
    print(f"  sat_id_dir: {sat_id_dir}")
    print(f"  rmin: {rmin}")
    print(f"  rmax: {rmax}")
    print(f"  bmin: {bmin}")
    print(f"  bmax: {bmax}")
    print(f"  overdensity: {overdensity}")
    print(f"  sat_mask: {sat_mask}")
    print(f"  only_sat_mask: {only_sat_mask}")
    print(f"  halo_mask: {halo_mask}")
    print(f"  rand_mask: {rand_mask}")
    print(f"  subsample: {subsample}")
    print(f"  ananke_name: {ananke_name}")
    print(f"  photo_sys: {photo_sys}")
    print(f"  cmd_magnames: {cmd_magnames}")
    print(f"  app_mag_lim_lo: {app_mag_lim_lo}")
    print(f"  app_mag_lim_hi: {app_mag_lim_hi}")
    print(f"  abs_mag_lim_lo: {abs_mag_lim_lo}")
    print(f"  abs_mag_lim_hi: {abs_mag_lim_hi}")
    print(f"  ananke_r_max: {ananke_r_max}")
    
    
    
     
    # Read in star particles from the selected fire2 snapshot
    #redshift = 0.50239229 #snap 385
    redshift = 0.
    
    part = ga.io.Read.read_snapshots(species='star',
                                 snapshot_value_kind='redshift',
                                 snapshot_values=redshift,
                                 simulation_directory=sim_dir,
                                 elements='all',
                                 assign_hosts=True,
                                 assign_hosts_rotation=True,
                                 assign_orbits=True)
    
    
    
    
    # array of positions in principal axes (in kpc)
    pos_pa = part['star'].prop('host.distance.principal')

    # array of velocities in principal axes (in km/s)
    vel_pa = part['star'].prop('host.velocity.principal')

    # array of star particle mass in solar masses
    mass = part['star']['mass']

    # array of decimal log stellar ages (in Gyr)
    log_age = np.log10(part['star'].prop('age') * 1e9)

    # array of star particle metallicities
    feh = part['star'].prop('metallicity.fe')

    abundances_list = ['helium', 'carbon', 'nitrogen', 'oxygen', 'neon', 'magnesium', 'silicon', 'sulfur', 'calcium']
    # dictionary of chemical abundance arrays (X/H)
    abundances = {'sulphur' if el == 'sulfur' else el: part['star'].prop('metallicity.' + el) for el in abundances_list}

    # alpha abundance (Mg/Fe)
    alpha = abundances['magnesium'] - feh
    
    
    
    
    # Preparing the star particles data to be used by ananke, masking only a sphere within 30 kpc
    #mask = np.linalg.norm(pos_pa, axis=1)<=30
    distance_mask = (np.linalg.norm(pos_pa, axis=1) >= rmin) & (np.linalg.norm(pos_pa, axis=1) <= rmax)
    
    # Blank satellite mask
    satellite_mask = np.ones(len(pos_pa), dtype=bool)
    
    # Store satellite indices (from Mia)
    all_unique_lmc_inds = np.loadtxt(sat_id_dir)
    all_inds = part['star'].prop('id')
    remaining_indices = np.setdiff1d(all_inds, np.array(all_unique_lmc_inds))
    id_indices = np.where(np.isin(all_inds, remaining_indices))[0]
    only_lmc_indices = np.where(~np.isin(all_inds, remaining_indices))[0]
    
    # Remove satellite
    if sat_mask:
        print("* Removing satellite with mask \n")
        
        # Update satellite mask 
        satellite_mask[only_lmc_indices] = False
        
   # Combine masks
    mask = distance_mask & satellite_mask
        
              

    p = {}
    p['pos3'] = pos_pa[mask]    # position in kpc (Nx3)
    p['vel3'] = vel_pa[mask]    # velocity in km/s (Nx3)
    p['mass'] = mass[mask]      # mass in solar masses
    p['age'] = log_age[mask]    # log age in Gyr
    p['feh'] = feh[mask]        # [Fe/H]
    for el, abun in abundances.items():  
        p[el] = abun[mask]      # other abundances as [X/H]

    p['alpha'] = alpha[mask]    # alpha abundance [Mg/Fe]

    p['parentid'] = np.where(mask)[0]           # indices of parent particles in snapshot
    p['dform'] = 0*p['mass']  # dummy variable for now
    
    
    
    
    # Now we can prepare the ananke surveyor. Default surveyor is set to simulate a Roman + HST photometric system.
    fsample = 0.005
    
    surveyor = an.Ananke(p, name=ananke_name, fsample=fsample, rSun0=0, rSun1=0, rSun2=0, 
                         photo_sys=photo_sys, cmd_magnames=cmd_magnames, abs_mag_lim_hi=abs_mag_lim_hi)
    
    survey = surveyor.run()
  