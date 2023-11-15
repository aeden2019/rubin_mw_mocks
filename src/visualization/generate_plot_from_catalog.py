"""
Script to generate mollweide plot from ananke mock catalog

"""

import numpy as np
import mollweide_plotting as pl
import nba
import vaex

if __name__ == "__main__":
    
    # Open the mock file
    df = vaex.open("ananke_mock.hdf5")

    # Combine the position and velocities
    pos_out = np.column_stack((df['px'].values, df['py'].values, df['pz'].values))
    vel_out = np.column_stack((df['vx'].values, df['vy'].values, df['vz'].values))
    
    # Apply kinematics 
    f = 1* (u.km/u.s).to(u.kpc/u.Gyr)
    kinematics1_ananke = nba.kinematics.Kinematics(pos_out, vel_out)
    pos_galactic_ananke = kinematics1_ananke.pos_cartesian_to_galactic()
    vel_galactic_ananke = kinematics1_ananke.vel_cartesian_to_galactic()
    
    # Create mollwiede plot
    figname = "mock_mollweide_plot.png"
    pl.mollweide_projection(pos_galactic_ananke[0]*180/np.pi, pos_galactic_ananke[1]*180/np.pi, 0, 0, 
                            sim_dir=sim_dir, bmin=bmin, bmax=bmax, nside=40, smooth=5, overdensity=overdensity, figname=figname)