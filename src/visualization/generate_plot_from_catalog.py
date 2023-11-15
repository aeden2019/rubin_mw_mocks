"""
Script to generate mollweide plot from ananke mock catalog

"""

import numpy as np
import mollweide_plotting as pl
import nba
import vaex
from astropy import units as u

import sys
sys.path.append('../mock')
import allvariables

if __name__ == "__main__":
  
    # Get the necessary variables from the parameter file
    params = allvariables.getparams()
    sim_dir = params[2]
    bmin = params[5]
    bmax = params[6]
    overdensity = params[7]
    
    # Open the mock file
    mock_file_path = "/home/jovyan/home/rubin_mw_mocks/src/mock/ananke_mock.hdf5" # UPDATE THIS TO BE INDEPENDENT
    print(f"Opening data from: {mock_file_path}")
    df = vaex.open(mock_file_path)

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
    
    # Create mollwiede plot
    print("Creating mollweide plot")
    figname = "mock_mollweide_plot.png"
    pl.mollweide_projection(pos_galactic_ananke[0]*180/np.pi, pos_galactic_ananke[1]*180/np.pi, 0, 0, 
                            sim_dir=sim_dir, bmin=bmin, bmax=bmax, nside=40, smooth=5, overdensity=overdensity, figname=figname)