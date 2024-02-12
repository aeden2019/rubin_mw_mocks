"""
Satellite removal functions from Arpit. 

Note, you want to use the new rockstar catalogs i.e rockstar_dm_new to do all of this.

Now the following fx call for whatever simulation directory you want and snapshot  will remove all the substructure with total mass between 1e8--3e11:
MW_dark_inds = cc.remove_dark_substructure(sim_dir, nsnap, mass_cuts=[1e8, 3e11])

You can plot them as:
pos_dark = part['dark'].prop('host.distance.principal')
plt.hist2d(pos_dark[MW_dark_inds,0], pos_dark[MW_dark_inds,1], bins=np.linspace(-300, 300, 300), norm=LogNorm(),)

"""



## Reading Rockstar dark matter particle assignments. 
def read_particle_assignments(sim_dir, nsnap, rockstar_directory='halo/rockstar_dm_new/'):
    """
    Reads dark matter particle assignments from Rockstar arranged by halo_ID: [particle IDs]
    dictionary structure.

    Args:
        sim_dir (str): The simulation directory containing the Rockstar particle assignment data.
        nsnap (int): The snapshot number to read the particle assignments for.
        rockstar_directory (str, optional): The subdirectory path within `sim_dir` where Rockstar data is stored.
                                            Default is 'halo/rockstar_dm_new/'.

    Returns:
        dict: A dictionary mapping halo_IDs to lists of particle IDs belonging to each halo.

    Note:
        This function reads the dark matter particle assignments computed by Rockstar for a specific snapshot.
        The particle assignments table should have been precomputed using a C++ pipeline and saved as a pickle file.
        The function will check if a corresponding pickle file exists in the given `sim_dir` before attempting to read
        the particle assignments from the original particle table. If the pickle file is not found, the function will
        generate it from the particle table and save it for faster access in future calls.

        To match particle IDs in the simulation snapshot, the particle ID for the dark matter species should be
        accessed as part['dark']['id'] from the Gizmo snapshot data.

    Examples:
        # Read particle assignments for snapshot 10 in the 'halo/rockstar_dm_new/' directory
        particle_assignments = read_particle_assignments(sim_dir="/path/to/sim_directory", nsnap=500)

    """
    import os, pickle
    
    file_path = f'{sim_dir}/{rockstar_directory}catalog_hdf5/dark_assignment_{nsnap:03d}.pickle'
    
    if os.path.isfile(file_path):
        print('Reading particle assignments from pickle file.')
        with open(file_path, 'rb') as f:
            return pickle.load(f)
        
    print('Generating pickle file from particle table')
    
    def process_line(line):
        halo_id, particle_ids = line.split(":")
        return {int(halo_id): list(map(int, particle_ids.split(",")[:-1]))}
    
    particle_file_dir = f'{sim_dir}/{rockstar_directory}catalog/particle_table:{nsnap:03d}'
    final_dict = {}
    
    with open(particle_file_dir, "r") as f:
        [final_dict.update(line_dict) for line_dict in map(process_line, f.read().strip().split("\n"))]
    
    print(f'Saving assignments @ {file_path}')
    
    try:
        with open(file_path, 'wb') as f:
            pickle.dump(final_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
        print('Files saved!')
    except Exception as e:
        print(e)
        print('Saving failed. Make sure you have writing permissions to the dir or change the save dir.')
        
    return final_dict





def remove_dark_substructure(sim_dir, nsnap,
                             part=None, hal_IDs=None,
                             mass_cuts=[1e9, 1e11],  # [lower_limit, upper_limit] on mass
                             return_sub=False, **kwargs):
    """
    Removes substructure by halo mass_limits using particle assignments from Rockstar.

    Note:
    - Particle assignments table should be precomputed using a C++ pipeline and provided.
    - If `hal_IDs` is not provided, the function will impose `mass_cuts` to determine the halos to remove.

    Args:
        sim_dir (str): The simulation directory where the required data is stored.
        nsnap (int): The snapshot index indicating which snapshot to process.
        part (dict, optional): Particle directory from the Gizmo snapshot. 
                               If not provided, it will be read from the simulation directory and snapshot index.
        hal_IDs (list, optional): A list of halo IDs. 
                                  If provided, the function will remove substructures associated with these halos 
                                  based on the mass limits. If not provided, the function will use `mass_cuts` 
                                  to determine the halo IDs to remove.
        mass_cuts (list, optional): A list containing two elements representing the mass limits 
                                    [lower_limit, upper_limit] on halo mass in solar masses.
        return_sub (bool, optional): If True, the function also returns the removed particle indices.
        **kwargs: Additional keyword arguments to be passed to underlying functions.

    Returns:
        np.array: An array containing the indices of DM particles after removing substructures 
                  by mass limits on halos.

    Examples:
        # Remove substructures with specified hal_IDs and get the remaining DM particle indices
        remaining_indices = remove_dark_substructure(sim_dir="/path/to/sim_directory", 
                                                     nsnap=500, 
                                                     hal_IDs=[1, 2, 3],
                                                     return_sub=False)

        # Remove substructures using mass_cuts and get the remaining DM particle indices and removed substructure indices
        remaining_indices, removed_indices = remove_dark_substructure(sim_dir="/path/to/sim_directory", 
                                                                      nsnap=500, 
                                                                      mass_cuts=[1e9, 1e11],
                                                                      return_sub=True)

    """
    import halo_analysis as halo
    import numpy as np
    
    if hal_IDs is None:
        hal = halo.io.IO.read_catalogs('index', nsnap, sim_dir, species=['dark'], **kwargs)
        mass_filt = np.where((hal['mass'] > mass_cuts[0]) & ((hal['mass'] < mass_cuts[1])))[0]
        hal_IDs = hal['id'][mass_filt]
    
    final_dict = read_particle_assignments(sim_dir, nsnap, **kwargs)
    dark_IDs = np.hstack([final_dict[key] for key in hal_IDs])
    
    if part is None:
        part = read(sim_dir, nsnap, **kwargs)
    
    remove_dark_inds = np.isin(part['dark']['id'], dark_IDs)
    
    if return_sub:
        return np.where(~remove_dark_inds)[0], np.where(remove_dark_inds)[0]##main halo, substrucutre
    
    return np.where(~remove_dark_inds)[0]



