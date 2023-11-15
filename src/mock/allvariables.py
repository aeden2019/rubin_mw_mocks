"""
Read variables from yaml parameter file

"""


import yaml
from argparse import ArgumentParser
import os



def getparams():
    # Parse the parameter file
    parser = ArgumentParser(description="Parameters file for make_mock.py")
    parser.add_argument(
            "--param", dest="paramFile", default="config.yaml",
            type=str, help="provide parameter file")
    args = parser.parse_args()
    
    # Pass the parameter file to readparams
    paramfile = args.paramFile
    params = readparams(paramfile)
    
    return params



def readparams(paramfile):
    
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construct the full path to the parameter file
    paramfile_path = os.path.join(script_dir, paramfile)
    
    # Open the param file path
    with open(paramfile_path) as f:
        d = yaml.safe_load(f)
    
    # Store each yaml entry
    snap = d["snap"]
    sim = d["sim"]
    sim_dir = d["sim_dir"]
    rmin = d["rmin"]
    rmax = d["rmax"]
    bmin = d["bmin"]
    bmax = d["bmax"]
    overdensity = d["overdensity"]
    sat_mask = d["sat_mask"]
    halo_mask = d["halo_mask"]
    rand_mask = d["rand_mask"]
    subsample = d["subsample"]
    
    # Check for correct type
    assert type(snap)==int, "snap must be an integer"
    assert type(sim)==str, "sim must be a string"
    assert type(sim_dir)==str, "sim_dir must be a string"
    assert type(rmin)==int, "rmin must be an integer"
    assert type(rmax)==int, "rmax must be an integer"
    assert type(bmin)==int, "bmin must be an integer"
    assert type(bmax)==int, "bmax must be an integer"
    assert type(overdensity) == bool, "overdensity must be a bool"
    assert type(sat_mask) == bool, "sat_mask must be a bool"
    assert type(halo_mask) == bool, "halo_mask must be a bool"
    assert type(rand_mask) == bool, "rand_mask must be a bool"
    assert type(subsample)==int, "subsample must be an integer"

    return [snap, sim, sim_dir, rmin, rmax, bmin, bmax, overdensity, sat_mask, halo_mask, rand_mask, subsample]