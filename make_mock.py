"""
Create ANANKE mock catalog from FIRE data

"""

# import gizmo_analysis as ga
# import sys
# import graphing 
# import ananke as an
# import numpy as np
# from astropy import units as u
# import nba
# import plotting as pl

import numpy as np
import args
import gizmo_analysis as ga
import allvariables
from argparse import ArgumentParser
from astropy import units as u
import nba
import ananke as an

# CREATE INDEPENDENT FILES FOR THESE TOOLS
# import sys
# sys.path.append('/home/jovyan/home/python_tools')

import mollweide_plotting as pl



if __name__ == "__main__":
    
    # Parse the parameter file
    parser = ArgumentParser(description="Parameters file for make_mock.py")
    parser.add_argument(
            "--param", dest="paramFile", default="config.yaml",
            type=str, help="provide parameter file")
    args = parser.parse_args()
    
    # Pass the parameter file to allvariables.readparams
    paramfile = args.paramFile
    params = allvariables.readparams(paramfile)
    
    # Declare all the parameters locally
    snap = params[0]
    sim = params[1]
    sim_dir = params[2]
    rmin = params[3]
    rmax = params[4]
    bmin = params[5]
    bmax = params[6]
    overdensity = params[7]
    sat_mask = params[8]
    halo_mask = params[9]
    rand_mask = params[10]
    subsample = params[11]
    
    
    # Create FIRE part
    part = ga.io.Read.read_snapshots(['star'], 'index', snap, sim_dir, assign_hosts=True)
    
    # Declare position and velocity
    pos = part['star'].prop('host.distance')
    vel = part['star'].prop('host.velocity')
    
    
    # Blank satellite mask
    satellite_mask = np.ones(len(pos), dtype=bool)
    
    # Remove satellite
    if sat_mask:
        # Store satellite indices (from Mia)
        path = '/home/jovyan/home/tracked/{}/all_LMC_inds.txt'.format(sim)
        all_unique_lmc_inds = np.loadtxt(path)
        all_inds = part['star'].prop('id')
        remaining_indices = np.setdiff1d(all_inds, np.array(all_unique_lmc_inds))
        id_indices = np.where(np.isin(all_inds, remaining_indices))[0]
        only_lmc_indices = np.where(~np.isin(all_inds, remaining_indices))[0]
        
        # Update satellite mask 
        satellite_mask[only_lmc_indices] = 0

        # Apply satellite mask
        pos = pos[satellite_mask]     
        vel = vel[satellite_mask]
         
            
    # Blank halo mask
    dist_cut1 = np.ones(len(pos), dtype=bool)
        
    # Removing outer halo
    if halo_mask:
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
    
    # Run the ananke process with kword args
    name='sim'
    ananke = an.Ananke(p, name, photo_sys='padova/LSST', cmd_magnames='rmag,gmag-rmag'
                                                , app_mag_lim_lo=17, app_mag_lim_hi=27.5, abs_mag_lim_lo=-7.0, abs_mag_lim_hi=10.0
                                                , color_lim_lo=-1000, color_lim_hi=1000, r_max=1000)
    ananke.run()
    
    # Make a survey using LSST (outputted as vaex data structure)
    survey = ananke._output
    df = survey._vaex
    
    # Combine the position and velocities
    pos_out = np.column_stack((df['px'].values, df['py'].values, df['pz'].values))
    vel_out = np.column_stack((df['vx'].values, df['vy'].values, df['vz'].values))
    
    # Apply kinematics 
    f = 1* (u.km/u.s).to(u.kpc/u.Gyr)
    kinematics1_ananke = nba.kinematics.Kinematics(pos_out, vel_out)
    pos_galactic_ananke = kinematics1_ananke.pos_cartesian_to_galactic()
    vel_galactic_ananke = kinematics1_ananke.vel_cartesian_to_galactic()
    
    # Create mollwiede plot
    # UPDATE pl TO BE A SEPERATE FILE
#     pl.mollweide_projection(pos_galactic_ananke[0]*180/np.pi, pos_galactic_ananke[1]*180/np.pi, 0, 0, 
#                             title="", bmin=bmin, bmax=bmax, nside=40, smooth=5, overdensity=overdensity, figname=figname)
    figname = "mock_mollweide_plot.png"
    pl.mollweide_projection(pos_galactic_ananke[0]*180/np.pi, pos_galactic_ananke[1]*180/np.pi, 0, 0, 
                            sim_dir=sim_dir, bmin=bmin, bmax=bmax, nside=40, smooth=5, overdensity=overdensity, figname=figname)
    
    
    
    