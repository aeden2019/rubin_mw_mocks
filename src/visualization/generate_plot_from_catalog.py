"""
Script to generate mollweide overdensity plot from ananke mock catalog

"""

import numpy as np
import mollweide_plotting as pl
import nba
import vaex
from astropy import units as u

import sys
sys.path.append('../mock')
import allvariables
import os

if __name__ == "__main__":
  
    # Get the parameters from the config file
    params = allvariables.getparams()

    # Define parameter names
    param_names = [
        "snap", "sim", "sim_dir", "sat_id_dir", "rmin", "rmax", "sat_mask", "ananke_name", 
        "photo_sys", "cmd_magnames", "app_mag_lim_lo", "app_mag_lim_hi", 
        "abs_mag_lim_lo", "abs_mag_lim_hi", "ananke_r_max", "fsample"
    ]

    # Initialize parameter dictionary
    params_dict = {}

    # Assign values to parameters and declare them locally
    for param_name, param_value in zip(param_names, params):
        params_dict[param_name] = param_value
        locals()[param_name] = param_value
    
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the full path to the data file
    ananke_file_path = os.path.join(script_dir, "..", "mock", f"survey.{ananke_name}.h5")
    print(f"\nOpening data from: {ananke_file_path}")
    df = vaex.open(ananke_file_path)
    
    # Combine the position and velocities
    print("Combining position and velocity")
    pos_out = np.column_stack((df['px'].values, df['py'].values, df['pz'].values))
    vel_out = np.column_stack((df['vx'].values, df['vy'].values, df['vz'].values))
    
    # Apply kinematics 
    print("Applying kimematics")
    f = 1* (u.km/u.s).to(u.kpc/u.Gyr)
    kinematics1_ananke = nba.kinematics.Kinematics(pos_out, vel_out)
    pos_galactic_ananke = kinematics1_ananke.pos_cartesian_to_galactic()
    vel_galactic_ananke = kinematics1_ananke.vel_cartesian_to_galactic()
    
    # Declare default params for overdensity
    bmin = -1
    bmax = 1
    overdensity = True
    
    # Create mollwiede overdensity plot
    print("Creating mollweide overdensity plot")
    figname = f"{ananke_name}_mollweide_overdensity_plot.png"
    pl.mollweide_projection(pos_galactic_ananke[0]*180/np.pi, pos_galactic_ananke[1]*180/np.pi, 0, 0, 
                            sim_dir=sim_dir, bmin=bmin, bmax=bmax, nside=40, smooth=5, overdensity=overdensity, figname=figname)
    