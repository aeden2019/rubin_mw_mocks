"""
Read variables from yaml parameter file

"""


import yaml
from argparse import ArgumentParser
import os


def getparamsnotebook(paramfile="config.yaml"):
    """
    Parses the config.yaml file for use in readparams(). 
    This version is specific for notebooks to avoid command line argument parsing. 
    
    Output:
    paramfile = parsed config.yaml file from getparams()
    
    """
    
    params = readparams(paramfile)
    
    return params


def getparams():
    """
    Parses the config.yaml file for use in readparams(). 
    
    Output:
    paramfile = parsed config.yaml file from getparams()
    
    """
    
    # Parse the parameter file
    parser = ArgumentParser(description="Parameters file for make_mock.py")
    parser.add_argument(
            "--param", dest="paramFile", default="config.yaml",
            type=str, help="provide parameter file")
    arguments = parser.parse_args()
    
    # Pass the parameter file to readparams
    paramfile = arguments.paramFile
    params = readparams(paramfile)
    
    return params



def readparams(paramfile):
    """
    Reads the output of getparams(), asserts that each parameter is the correct type, and returns an array of parameters. 
     
    Input: 
    paramfile = parsed config.yaml file from getparams()
    
    Output:
    array of parameters in correct types
    
    """
    
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
    sat_id_dir = d["sat_id_dir"]
    rmin = d["rmin"]
    rmax = d["rmax"]
    sat_mask = d["sat_mask"]
    ananke_name = d["ananke_name"]
    photo_sys = d["photo_sys"]
    cmd_magnames = d["cmd_magnames"]
    app_mag_lim_lo = d["app_mag_lim_lo"]
    app_mag_lim_hi = d["app_mag_lim_hi"]
    abs_mag_lim_lo = d["abs_mag_lim_lo"]
    abs_mag_lim_hi = d["abs_mag_lim_hi"]
    ananke_r_max = d["ananke_r_max"]
    fsample = d["fsample"]
    
    # Check for correct type
    assert isinstance(snap, int), f"snap must be of type int, but got {type(snap)}"
    assert isinstance(sim, str), f"sim must be of type str, but got {type(sim)}"
    assert isinstance(sim_dir, str), f"sim_dir must be of type str, but got {type(sim_dir)}"
    assert isinstance(sat_id_dir, str), f"sat_id_dir must be of type str, but got {type(sat_id_dir)}"
    assert isinstance(rmin, int), f"rmin must be of type int, but got {type(rmin)}"
    assert isinstance(rmax, int), f"rmax must be of type int, but got {type(rmax)}"
    assert isinstance(sat_mask, bool), f"sat_mask must be of type bool, but got {type(sat_mask)}"
    assert isinstance(ananke_name, str), f"ananke_name must be of type str, but got {type(ananke_name)}"
    assert isinstance(photo_sys, str), f"photo_sys must be of type str, but got {type(photo_sys)}"
    assert isinstance(cmd_magnames, str), f"cmd_magnames must be of type str, but got {type(cmd_magnames)}"
    assert isinstance(app_mag_lim_lo, float), f"app_mag_lim_lo must be of type float, but got {type(app_mag_lim_lo)}"
    assert isinstance(app_mag_lim_hi, float), f"app_mag_lim_hi must be of type float, but got {type(app_mag_lim_hi)}"
    assert isinstance(abs_mag_lim_lo, float), f"abs_mag_lim_lo must be of type float, but got {type(abs_mag_lim_lo)}"
    assert isinstance(abs_mag_lim_hi, float), f"abs_mag_lim_hi must be of type float, but got {type(abs_mag_lim_hi)}"
    assert isinstance(ananke_r_max, int), f"ananke_r_max must be of type int, but got {type(ananke_r_max)}"
    assert isinstance(fsample, float), f"fsample must be of type float, but got {type(fsample)}"



    return [snap, sim, sim_dir, sat_id_dir, rmin, rmax, sat_mask, ananke_name,
            photo_sys, cmd_magnames, app_mag_lim_lo, app_mag_lim_hi,
            abs_mag_lim_lo, abs_mag_lim_hi, ananke_r_max, fsample]