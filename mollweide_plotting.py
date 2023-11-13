"""
Mollweide plotting routine, addapted from plotting.py written by Nico Garavito-Camargo (jngaravitoc)

"""


import numpy as np
import matplotlib.pyplot as plt
import pynbody
import sys
from pynbody import filt
from astropy import units as u
# from astropy.coordinates import Latitude
# from astropy.coordinates import Longitude
import nba
import healpy as hp
from healpy.newvisufunc import projview, newprojplot



def mollweide_projection(l, b, l2, b2, sim_dir, bmin, bmax, nside, smooth, q=[0], **kwargs):

    """
    Makes mollweide plot using healpix
    Parameters:
    ----------- 
    l : numpy.array in degrees 
    b : numpy.array in degrees [-90, 90]
    """
 
    times = '{}snapshot_times.txt'.format(sim_dir)
    
    mwlmc_indices = hp.ang2pix(nside,  (90-b)*np.pi/180., l*np.pi/180.)
    npix = hp.nside2npix(nside)
 
    idx, counts = np.unique(mwlmc_indices, return_counts=True)
    degsq = hp.nside2pixarea(nside, degrees=True)
    # filling the full-sky map
    hpx_map = np.zeros(npix, dtype=float)
    if q[0] != 0 :    
        counts = np.zeros_like(idx, dtype=float)
        k=0
        for i in idx:
            pix_ids = np.where(mwlmc_indices==i)[0]
            counts[k] = np.mean(q[pix_ids])
            k+=1
        hpx_map[idx] = counts
    else :
       hpx_map[idx] = counts/degsq
   
    
    # Check for overdensity
    map_smooth = hp.smoothing(hpx_map, fwhm=smooth*np.pi/180)
    overdensity = False
    if 'overdensity' in kwargs.keys():
        overdensity = kwargs['overdensity']

        # Overdensity calculation
        if overdensity == True:
            map_smooth = (map_smooth / np.mean(map_smooth)) - 1
    
    if 'cmap' in kwargs.keys():
        cmap = kwargs['cmap']
    else:
        cmap='viridis'
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    plt.close()
    
    projview(
      map_smooth,
      coord=["G"],
      graticule=True,
      graticule_labels=True,
      unit=" ",
      xlabel="Longitude (l) ",
      ylabel="Latitude (b)",
      cb_orientation="horizontal",
      min=bmin,
      max=bmax,
      latitude_grid_spacing=45,
      projection_type="mollweide",
      cmap=cmap,
      fontsize={
              "xlabel": 25,
              "ylabel": 25,
              "xtick_label": 20,
              "ytick_label": 20,
              "cbar_label": 20,
              "cbar_tick_label": 20,
              }
      )
    

    # CHANGED COLOR TO RED FROM YELLOW
    newprojplot(theta=np.radians(90-(b2)), phi=np.radians(l2), marker="o", color="red", markersize=5, lw=0, mfc='none')
    if 'l3' in kwargs.keys():
        l3 = kwargs['l3']
        b3 = kwargs['b3']
        newprojplot(theta=np.radians(90-(b3)), phi=np.radians(l3), marker="o", color="red", markersize=5, lw=0)
    elif 'l4' in kwargs.keys():
        l4 = kwargs['l4']
        b4 = kwargs['b4']
        newprojplot(theta=np.radians(90-(b4)), phi=np.radians(l4), marker="*", color="r", markersize=8, lw=0)

    #newprojplot(theta=np.radians(90-(b2[0])), phi=np.radians(l2[0]-120), marker="*", color="r", markersize=5 )
    #newprojplot(theta=np.radians(90-(b2[1])), phi=np.radians(l2[1]-120), marker="*", color="w", markersize=2 )
    
    # TESTING - MOVE COLORBAR DOWN AND ADD LABEL 
    cax = plt.gcf().get_axes()[-1]
    cax.set_position([cax.get_position().x0, cax.get_position().y0 - 0.1, cax.get_position().width, cax.get_position().height])
    
    # TESTING - UPDATE LABEL BASED ON SIMULATION TYPE
    if 'velocity' in kwargs.keys():
        velocity = kwargs['velocity']
        if velocity:
            if overdensity:
                cax.set_xlabel(r"$\Delta|v|$ / $|v|$", fontsize=20)
            else:
                cax.set_xlabel(r"|$v$| [km/s]", fontsize=20)
        else: 
            if overdensity: 
                cax.set_xlabel(r"$\Delta\rho$ / $\rho$", fontsize=20)
            else: 
                cax.set_xlabel(r"$\rho$", fontsize=20)
    else: 
        if overdensity: 
            cax.set_xlabel(r"$\Delta\rho$ / $\rho$", fontsize=20)
        else: 
            cax.set_xlabel(r"$\rho$", fontsize=20)
       
    
    if 'figname' in kwargs.keys():
        print("* Saving figure in ", kwargs['figname'])
        plt.savefig(kwargs['figname'], bbox_inches='tight')
        plt.close()
        
    plt.close()
    #return 0