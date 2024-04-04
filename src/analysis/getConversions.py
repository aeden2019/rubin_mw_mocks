"""
Functions that return commonly used conversions. 


"""

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



