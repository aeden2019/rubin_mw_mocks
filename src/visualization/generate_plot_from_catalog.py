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
import os

if __name__ == "__main__":
  
    # Get the necessary variables from the parameter file
    params = allvariables.getparams()
    sim_dir = params[2]
    bmin = params[5]
    bmax = params[6]
    overdensity = params[7]
    
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the full path to the data file
    ananke_file_path = os.path.join(script_dir, "..", "mock", "survey.sim.h5")
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
    
    # Create mollwiede plot
    print("Creating mollweide plot")
    figname = "ananke_mollweide_plot.png"
    pl.mollweide_projection(pos_galactic_ananke[0]*180/np.pi, pos_galactic_ananke[1]*180/np.pi, 0, 0, 
                            sim_dir=sim_dir, bmin=bmin, bmax=bmax, nside=40, smooth=5, overdensity=overdensity, figname=figname)
    