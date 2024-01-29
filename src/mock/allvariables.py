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
    sat_id_dir = d["sat_id_dir"]
    rmin = d["rmin"]
    rmax = d["rmax"]
    bmin = d["bmin"]
    bmax = d["bmax"]
    overdensity = d["overdensity"]
    sat_mask = d["sat_mask"]
    only_sat_mask = d["only_sat_mask"]
    halo_mask = d["halo_mask"]
    rand_mask = d["rand_mask"]
    subsample = d["subsample"]
    ananke_name = ["ananke_name"]
    photo_sys = ["photo_sys"]
    cmd_magnames = ["cmd_magnames"]
    app_mag_lim_lo = ["app_mag_lim_lo"]
    app_mag_lim_hi = ["app_mag_lim_hi"]
    abs_mag_lim_lo = ["abs_mag_lim_lo"]
    abs_mag_lim_hi = ["abs_mag_lim_hi"]
    ananke_r_max = ["ananke_r_max"]
    
    # Check for correct type
    assert type(snap)==int, "snap must be an integer"
    assert type(sim)==str, "sim must be a string"
    assert type(sim_dir)==str, "sim_dir must be a string"
    assert type(sat_id_dir)==str, "sat_id_dir must be a string"
    assert type(rmin)==int, "rmin must be an integer"
    assert type(rmax)==int, "rmax must be an integer"
    assert type(bmin)==int, "bmin must be an integer"
    assert type(bmax)==int, "bmax must be an integer"
    assert type(overdensity) == bool, "overdensity must be a bool"
    assert type(sat_mask) == bool, "sat_mask must be a bool"
    assert type(only_sat_mask) == bool, "only_sat_mask must be a bool"
    assert type(halo_mask) == bool, "halo_mask must be a bool"
    assert type(rand_mask) == bool, "rand_mask must be a bool"
    assert type(subsample)==int, "subsample must be an integer"
    assert type(ananke_name)==str, "ananke_name must be a string"
    assert type(photo_sys)==str, "photo_sys must be a string"
    assert type(cmd_magnames)==str, "cmd_magnames must be a string"
    assert type(app_mag_lim_lo)==float, "app_mag_lim_lo must be a float"
    assert type(app_mag_lim_hi)==float, "app_mag_lim_hi must be a float"
    assert type(abs_mag_lim_lo)==float, "abs_mag_lim_lo must be a float"
    assert type(abs_mag_lim_hi)==float, "abs_mag_lim_hi must be a float"
    assert type(ananke_r_max)==int, "ananke_r_max must be an integer"
    
    # Check that sat_mask and only_sat_mask are not both True
    assert not (sat_mask and only_sat_mask), "sat_mask and only_sat_mask cannot both be True"

    return [snap, sim, sim_dir, sat_id_dir, rmin, rmax, bmin, bmax, overdensity,
            sat_mask, only_sat_mask, halo_mask, rand_mask, subsample, ananke_name,
            ananke_name, photo_sys, cmd_magnames, app_mag_lim_lo, app_mag_lim_hi,
            abs_mag_lim_lo, abs_mag_lim_hi, ananke_r_max]