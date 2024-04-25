"""
Functions that return commonly used conversions. 


"""
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np




def abs2app(M, d, A=0.0):
    """
    Converts absolute magnitude to apparent magnitude.
    
    Parameters:
    -----------
    M : float
        Absolute magnitude (mag).
    d : float
        Distance (pc). 
    A : float
        Extinction (mag/pc).
    
    
    Returns:
    --------
    m : float
        Apparent magnitude (mag).  
    """
    
    m = 5 * np.log10(d) - 5 + M + A*d
    
    return m




def app2abs(m, d, A=0.0):
    """
    Converts apparent magnitude to absolute magnitude.
    
    Parameters:
    -----------
    m : float
        Apparent magnitude (mag).
    d : float
        Distance (pc). 
    A : float
        Extinction (mag/pc).
    
    
    Returns:
    --------
    M : float
        Absolute magnitude (mag).  
    """
    
    M = m - 5 * np.log10(d) + 5 - A*d
    
    return M




def xyz_to_galactic(x, y, z):
    """
    Converts cartesian coordinates to galactic coordinates in a galactocentric model. 
    
    Parameters:
    -----------
    x : float
        x position (kpc).
    y : float
        y position (kpc).
    z : float
        z position (kpc).
    
    
    Returns:
    --------
    l : float
        Galactic longitude (degrees).
    b : float 
        Galactic latitude (degrees).
    """
    
    # Create SkyCoord object from XYZ coordinates
    xyz = SkyCoord(x=x*u.kpc, y=y*u.kpc, z=z*u.kpc, frame='galactocentric', representation_type='cartesian')
    
    # Transform XYZ to Galactic coordinates (l, b, distance)
    galactic = xyz.transform_to('galactic')
    
    return galactic.l.deg, galactic.b.deg



