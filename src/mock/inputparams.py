"""
Script to generate plots for the distribution of stellar 
ages and metallicity of the particles fed into py-ananke

"""

import numpy as np
import matplotlib.pyplot as plt

def plot_stellar_ages(ages):

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
    
    # Create a histogram of stellar ages
    plt.figure(figsize=(8, 6))
    plt.hist(feh, bins=500)

    plt.xlabel("metalicities [Fe/H]")
    plt.ylabel("frequency")

    # Show the plot
    plt.show()
    
    # Save the plot
    plt.savefig("ananke_input_metalicities.png")
    