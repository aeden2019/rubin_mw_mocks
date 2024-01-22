"""
Script to generate plots for the distribution of stellar 
ages and metallicity of the particles fed into py-ananke

"""

import numpy as np
import matplotlib.pyplot as plt

def plot_stellar_ages(ages):
    """
    Creates a histogram of stellar ages and saves them to ananke_input_ages.png
    
    Input:
    ages = array of stellar ages in [log Gyr]
    
    """

    # Create a histogram of stellar ages
    plt.figure(figsize=(8, 6))
    plt.hist(ages, bins=500)

    plt.xlabel("ages [log age in Gyr]")
    plt.ylabel("frequency")

    # Show the plot
    plt.show()
    
    # Save the plot
    plt.savefig("ananke_input_ages.png")
    

def plot_metalicity(feh):
    """
    Creates a histogram of metalicity and saves them to ananke_input_metalicities.png
    
    Input:
    feh = array of metalicities in [Fe/H]
    
    """
    
    # Create a histogram of stellar ages
    plt.figure(figsize=(8, 6))
    plt.hist(feh, bins=500)

    plt.xlabel("metalicities [Fe/H]")
    plt.ylabel("frequency")

    # Show the plot
    plt.show()
    
    # Save the plot
    plt.savefig("ananke_input_metalicities.png")
    