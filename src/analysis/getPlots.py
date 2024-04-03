"""
Functions for creating regularly used plots, such as color magnitude diagrams. 

"""

import matplotlib
from matplotlib import pyplot as plt



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
    
    
    