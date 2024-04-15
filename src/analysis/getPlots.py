"""
Functions for creating regularly used plots, such as color magnitude diagrams. 

"""

import matplotlib
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
import numpy as np




def colormagplot(xaxis, yaxis):
    """
    Create color magnitude diagram.
    
    Parameters:
    -----------
    xaxis : np.array
        Array with g-r data.
    yaxis : np.array
        Array with M_g data.
    
    
    Returns:
    --------
    
    """
    
    # Define hex cmap
    hex_cmap = 'cividis'
    hex_cmap = plt.get_cmap(hex_cmap)

    # Create the subplot
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))

    # Create the hexbin, with 241 log bins and the hex_cmap colormap
    hb = axs.hexbin(xaxis, yaxis, gridsize=(241), bins='log', cmap=hex_cmap)

    # Setting labels
    axs.set_xlabel('$g - r$')
    axs.set_ylabel('$M_g$')
    axs.set_aspect(1.0/axs.get_data_ratio())  # Set aspect ratio
    axs.invert_yaxis()

    # Colorbar
    cbar = fig.colorbar(hb, ax=axs)
    cbar.set_label('Number of stars per hexcell')

    plt.tight_layout()
    plt.show()
    
    
    
    
def densityParams(radial, bins, rmin, rmax):
    """
    Calculates parameters used to plot the density profile.

    Parameters:
    -----------
    radial : numpy.ndarray
        Array of radial values.
    bins : integer
        Number of bins.
    rmin : integer
        Minimum radius in kpc. 
    rmax : integer
        Maximum radius in kpc. 

    Returns:
    --------
    numpy.ndarray, numpy.ndarray
        edges: radial locations of bin edges
        density: density values of each bin
    """
    
    # Set up the counts and edges based on limits
    counts, edges = np.histogram(radial, bins=bins, range=(rmin, rmax))
    
    # Calculate the volume and density
    volume = 4.0/3.0 * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
    density = counts / volume
    
    return edges, density




def interpolate_density(x, y):
    """
    Interpolates density values using linear interpolation.

    Parameters:
    -----------
    x : numpy.ndarray
        Array of x-values.
    y : numpy.ndarray
        Array of y-values.

    Returns:
    --------
    numpy.ndarray, numpy.ndarray
        Interpolated x-values and interpolated y-values.
    """
    # Remove zero values before interpolation
    nonzero_indices = y > 0
    x_nonzero = x[nonzero_indices]
    y_nonzero = y[nonzero_indices]
    
    # Perform linear interpolation
    linear_interp = interp1d(x_nonzero, y_nonzero, kind='linear', fill_value="extrapolate")
    
    # Interpolate values for all x
    y_interp = linear_interp(x)
    
    return x, y_interp




def interpolate_velocity(x, y):
    """
    Interpolates velocity values using linear interpolation.

    Parameters:
    -----------
    x : numpy.ndarray
        Array of x-values.
    y : numpy.ndarray
        Array of y-values.

    Returns:
    --------
    numpy.ndarray, numpy.ndarray
        Interpolated x-values and interpolated y-values.
    """
    # Remove zero and non-finite values before interpolation
    valid_indices = np.logical_and(y != 0, np.isfinite(y))
    x_valid = x[valid_indices]
    y_valid = y[valid_indices]
    
    # Perform linear interpolation
    linear_interp = interp1d(x_valid, y_valid, kind='linear', fill_value="extrapolate")
    
    # Interpolate values for all x
    y_interp = linear_interp(x)
    
    return x, y_interp




def normalize_data(data_interp):
    """
    Normalize the interpolated data by the value at the first index.
    
    Parameters:
    -----------
    data_interp : numpy.ndarray
        Interpolated data array.
    
    Returns:
    --------
    numpy.ndarray
        Normalized data array.
    """
    
    inner_value = data_interp[0]
    
    return data_interp / inner_value



