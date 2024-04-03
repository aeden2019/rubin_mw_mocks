"""
Functions that return the errors in parallax and proper motion measurements. 

The equations can be found in section 3.2 of Ivezic 2019, https://arxiv.org/pdf/0805.2366.pdf. 

"""

import getPhotometric
import numpy as np




def getSigmaAstro(colorband, theta=700.0, tvis=30.0, X=1.2)
    """
    Returns the total astrometric error.
    
    Parameters:
    -----------
    colorband : string
        Indicates the color band to use among 'ugrizy'.
    theta : float
        Random astrometric error of a visit. 
    tvis : float
        Exposure time in seconds. 
    X : float
        Airmass. 
    
    
    
    Returns:
    --------
    float 
        Total astrometric error.
    
    """
    
    # Get random and systematic errors
    sigmaRand = getRandAstro(colorband, theta, tvis, X)
    sigmaSyst = getSystAstro()
    
    sigma = np.sqrt(sigmaRand**2 + sigmaSyst**2)
    
    return sigma



def getRandAstro(colorband, theta=700.0, tvis=30.0, X=1.2):
    """
    Returns the random astrometric error.
    
    Parameters:
    -----------
    colorband : string
        Indicates the color band to use among 'ugrizy'.
    theta : float
        Random astrometric error of a visit. 
    tvis : float
        Exposure time in seconds. 
    X : float
        Airmass. 
    
    
    
    Returns:
    --------
    float 
        Random astrometric error.
    
    """
    
    Cm = getPhotometric.getCm(colorband, tvis)
    
    snr = getPhotometric.getm5(colorband, Cm, tvis, X)
    
    sigma = theta / snr
    
    return sigma




def getSystAstro():
    """
    Returns the systematic astrometric error.
    
    Parameters:
    -----------
    
    Returns:
    --------
    float 
        Systematic astrometric error.
    
    """
    
    sigma = 10.0
    
    return sigma



