"""
Create ANANKE mock catalog from FIRE data

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
    
    # Declare all the parameters locally
    snap = params[0]
    sim = params[1]
    sim_dir = params[2]
    sat_id_dir = params[3]
    rmin = params[4]
    rmax = params[5]
    bmin = params[6]
    bmax = params[7]
    overdensity = params[8]
    sat_mask = params[9]
    only_sat_mask = params[10]
    halo_mask = params[11]
    rand_mask = params[12]
    subsample = params[13]
    ananke_name = params[14]
    photo_sys = params[15]
    cmd_magnames = params[16]
    app_mag_lim_lo = params[17]
    app_mag_lim_hi = params[18]
    abs_mag_lim_lo = params[19]
    abs_mag_lim_hi = params[20]
    ananke_r_max = params[21]
    
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
     
    
    # Create FIRE part
    part = ga.io.Read.read_snapshots(['star'], 'index', snap, sim_dir, assign_hosts=True)
    
    # Declare position and velocity
    pos = part['star'].prop('host.distance')
    vel = part['star'].prop('host.velocity')
    
    
    # Blank satellite mask
    satellite_mask = np.ones(len(pos), dtype=bool)
    
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
        satellite_mask[only_lmc_indices] = 0

        # Apply satellite mask
        pos = pos[satellite_mask]     
        vel = vel[satellite_mask]
        
    # Include only the satellite mask
    if only_sat_mask:
        print("* Removing non-satellite objects with mask \n")
        
        # Update satellite mask 
        satellite_mask[id_indices] = 0

        # Apply satellite mask
        pos = pos[satellite_mask]     
        vel = vel[satellite_mask]
         
            
    # Blank halo mask
    dist_cut1 = np.ones(len(pos), dtype=bool)
        
    # Removing outer halo
    if halo_mask:
        print("* Removing halo with mask \n")
        # Create outer halo distance mask
        dist = np.sqrt(np.sum(pos**2, axis=1))
        dist_cut1 = np.where((dist > rmin) & (dist< rmax)) 
    
        # Apply outer halo distance mask
        pos = pos[dist_cut1]
        vel = vel[dist_cut1]
        
    
    # Blank random mask   
    rand_cut1 = np.ones(len(pos), dtype=bool)
   
    # Apply random subsampling
    if rand_mask:
        print(f"* Applying subsample of {subsample} particles to the dataset of {len(pos)} particles \n")
        # Create random cut
        rand_cut1 = np.random.randint(0, len(pos), subsample)
        pos = pos[rand_cut1]
        vel = vel[rand_cut1]
     
    #Create p dictionary to store particle data
    p = {}
    p['pos3'] = pos       # position in kpc
    p['vel3'] = vel       # velocity in km/s
    p['mass'] = part['star']['mass'][satellite_mask][dist_cut1][rand_cut1]                     # mass in solar masses
    p['age'] = part['star'].prop('age')[satellite_mask][dist_cut1][rand_cut1]                  # log age in Gyr
    p['feh'] = part['star'].prop('metallicity.fe')[satellite_mask][dist_cut1][rand_cut1]       # [Fe/H]
    p['helium'] = part['star'].prop('metallicity.he')[satellite_mask][dist_cut1][rand_cut1]    # [He/H]
    p['carbon'] = part['star'].prop('metallicity.c')[satellite_mask][dist_cut1][rand_cut1]     # [C/H]
    p['nitrogen'] = part['star'].prop('metallicity.n')[satellite_mask][dist_cut1][rand_cut1]   # [N/H]
    p['neon'] = part['star'].prop('metallicity.ne')[satellite_mask][dist_cut1][rand_cut1]      # [Ne/H]
    p['magnesium'] = part['star'].prop('metallicity.mg')[satellite_mask][dist_cut1][rand_cut1] # [Mg/H]
    p['silicon'] = part['star'].prop('metallicity.si')[satellite_mask][dist_cut1][rand_cut1]   # [Si/H]
    p['sulphur'] = part['star'].prop('metallicity.s')[satellite_mask][dist_cut1][rand_cut1]    # [S/H]
    p['calcium'] = part['star'].prop('metallicity.ca')[satellite_mask][dist_cut1][rand_cut1]   # [Ca/H]
    p['oxygen'] = part['star'].prop('metallicity.o')[satellite_mask][dist_cut1][rand_cut1]     # [O/H]
    p['alpha'] = part['star'].prop('metallicity.mg - metallicity.fe')[satellite_mask][dist_cut1][rand_cut1] # [Mg/Fe]
    p['parentid'] = part['star']['id'][satellite_mask][dist_cut1][rand_cut1] # indices of parent particles in snapshot
    p['dform'] = np.zeros(part['star']['position'].shape[0], dtype='float32')[satellite_mask][dist_cut1][rand_cut1] # dummy var
   
    # Assert that the lengths are the same
    assert len(p['pos3']) == len(p['vel3']), "Pos has different length than vel when creating p dictionary"
    assert len(p['pos3']) == len(p['mass']), "Pos has different length than mass when creating p dictionary"
    
    # Plotting for testing
    print("Plotting stellar ages and metalicity of ananke input")
    plot_ananke_inputs.plot_stellar_ages(p['age'])
    plot_ananke_inputs.plot_metalicity(p['feh'])
    
    fsample = 0.05
    
    # Initialize the ananke process with kword args
    ananke = an.Ananke(p, fsample=fsample , name=ananke_name, photo_sys=photo_sys, cmd_magnames=cmd_magnames, 
                       app_mag_lim_lo=app_mag_lim_lo, app_mag_lim_hi=app_mag_lim_hi,
                       abs_mag_lim_lo=abs_mag_lim_lo, abs_mag_lim_hi=abs_mag_lim_hi,
                       r_max=ananke_r_max)
    
    # Run ananke
    ananke.run()
   