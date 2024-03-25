"""
Functions that return the photometric error. 

The equations can be found in section 3.2 of Ivezic 2019, https://arxiv.org/pdf/0805.2366.pdf. 

"""

import numpy as np




def getSigmaVisit():
    """
    Returns the error of a single visit.
    
    Parameters:
    -----------
    colorband : string
        Indicates the color band to use among 'ugrizy'.
    mag : float
        Magnitude. 
    
    
    Returns:
    --------
    float
        Expected photometric error in magnitudes for a single visit. 
    """
    
    # Get the systematic and random errors
    sigma_sys = getSigmaSys()
    sigma_rand = getSigmaRand(colorband, mag)
    
    return np.sqrt(sigma_sys**2 + sigma_rand**2)




def getSigmaRand(colorband, mag, X=1.2, tvis=30):    
    """
    Returns the random photometric error.
    
    Parameters:
    -----------
    colorband : string
        Indicates the color band to use among 'ugrizy'.
    mag : float
        Magnitude. 
    tvis : float
        Exposure time in seconds. 
    X : float
        Airmass. 
    
    
    Returns:
    --------
    float
        Random photometric error.
    
    """
    
    # Get the parameter dictionary
    params = getEqnParameters(colorband)
    
    # Extract necessary parameters from the dictionary
    msky = params.get('msky')
    theta_eff = params.get('theta_eff')
    l = params.get('lambda')
    km = params.get('km')
    Cm = params.get('Cm')
    
    # Calculate the 5 sigma depth for point sources, eq. 6
    m_5 = Cm + 0.50*(msky - 21) + 2.5 * np.log(0.7 / theta_eff) + 1.25*np.log(tvis / 30) - km*(X - 1)
    
    # Calculate the x parameter
    x = 10**(0.4*(mag-m_5))
    
    # Calculate the random photmetric error for point sources, eq. 5
    sigma = np.sqrt((0.04 - l)*x + l*(x**2))
    
    return sigma




def getSigmaSys():
    """
    Returns the systematic photometric error.
    
    Parameters:
    -----------
    
    
    Returns:
    --------
    float
        Systematic photometric error.
        
    """
    
    # The sigma value must be below this value (MAKE THIS DYNAMIC)
    sigma = 0.005 
    
    return sigma
    
    

def getEqnParameters(colorband):
    """
    Returns the parameters in table 2 of Ivezic 2019.
    
    Parameters:
    -----------
    colorband : string
        Indicates the color band to use among 'ugrizy'.
    
    Returns:
    --------
    dictionary 
        'msky':        expected median zenith sky brightness at Cerro Pachon
        'theta':       expected delivered median zenith seeing
        'theta_eff':   effective zenith seeing  used for m5 computation
        'lambda':      band-dependent parameter from eq. 5
        'km':          adopted atmospheric extinction 
        'Cm':          band-dependent parameter from eq. 6
        'm5':          typical 5 sigma depth for point sources at zenith
        'delta_Cm':    loss of depth due to instrumental noise
        'delta_CM_2':  additive correction to Cm when exposure time is double to 60s
        'delta_m5':    loss depth at airmass of X = 1.2
        
    """
    
    # Check if colorband is valid
    valid_color_bands = ['u', 'g', 'r', 'i', 'z', 'y']
    if colorband not in valid_color_bands:
        raise ValueError("Invalid colorband. Please provide one of 'ugrizy'.")
    
    # Create a dictionary to store all the parameters
    data = {
        'u': {
            'msky': 22.99,
            'theta': 0.81,
            'theta_eff': 0.92,
            'lambda': 0.038,
            'km': 0.491,
            'Cm': 23.09,
            'm5': 23.78,
            'delta_Cm': 0.62,
            'delta_CM_2': 0.23,
            'delta_m5': 0.21
        },        
        'g': {
            'msky': 22.26,
            'theta': 0.77,
            'theta_eff': 0.87,
            'lambda': 0.039,
            'km': 0.213,
            'Cm': 24.42,
            'm5': 24.81,
            'delta_Cm': 0.18,
            'delta_CM_2': 0.08,
            'delta_m5': 0.16
        },
        'r': {
            'msky': 21.20,
            'theta': 0.73,
            'theta_eff': 0.83,
            'lambda': 0.039,
            'km': 0.126,
            'Cm': 24.44,
            'm5': 24.35,
            'delta_Cm': 0.10,
            'delta_CM_2': 0.05,
            'delta_m5': 0.14
        },
        'i': {
            'msky': 20.48,
            'theta': 0.71,
            'theta_eff': 0.80,
            'lambda': 0.039,
            'km': 0.096,
            'Cm': 24.32,
            'm5': 23.92,
            'delta_Cm': 0.07,
            'delta_CM_2': 0.03,
            'delta_m5': 0.13
        },
        'z': {
            'msky': 19.60,
            'theta': 0.69,
            'theta_eff': 0.78,
            'lambda': 0.039,
            'km': 0.069,
            'Cm': 24.16,
            'm5': 23.34,
            'delta_Cm': 0.05,
            'delta_CM_2': 0.02,
            'delta_m5': 0.13
        },        
        'y': {
            'msky': 18.61,
            'theta': 0.68,
            'theta_eff': 0.76,
            'lambda': 0.039,
            'km': 0.170,
            'Cm': 23.73,
            'm5': 22.45,
            'delta_Cm': 0.04,
            'delta_CM_2': 0.02,
            'delta_m5': 0.14
        }
    }
    
    return data.get(colorband)
    

