import os
import re
import h5py
import numpy as np
from glob import glob

SioverSil=0.163

def extract_rmax_zmax(name):
    """
    Extract rmax and zmax values (in kpc) from a depletion folder name.
    Example input: 'DEPLETION_rmax4.0kpc_zmax1.0kpc'
    Returns: (4.0, 1.0)
    """
    r_match = re.search(r"rmax([\d.]+)kpc", name)
    z_match = re.search(r"zmax([\d.]+)kpc", name)
    return float(r_match.group(1)), float(z_match.group(1))

def get_sim_name(path):
    """
    Extracts the name of the simulation from its full path.
    E.g., '/data/sim1' -> 'sim1'
    """
    return os.path.basename(os.path.abspath(path))

def load_depletion_file(filepath):
    """
    Loads a depletion file into a 2D NumPy array.
    Assumes the file is plain text with space-separated values.
    """
    return np.loadtxt(filepath)

def process_simulation(sim_path, h5file, output_number):
    """
    Processes one simulation directory:
    - Finds all depletion folders inside it.
    - Parses rmax/zmax from folder names.
    - Extracts only the specified depletion_output_XXXXX.txt file.
    - Stores the output as a dataset inside an HDF5 group.
    """
    sim_name = get_sim_name(sim_path)
    sim_group = h5file.require_group(sim_name)

    # Iterate through folders like DEPLETION_rmax..._zmax...
    for depletion_folder in os.listdir(sim_path):
        full_depl_path = os.path.join(sim_path, depletion_folder)
        if not os.path.isdir(full_depl_path) or not depletion_folder.startswith("DEPLETION_"):
            continue  # Skip non-folders or irrelevant names

        # Parse rmax and zmax from the folder name
        try:
            rmax, zmax = extract_rmax_zmax(depletion_folder)
        except Exception as e:
            print(f"Skipping {depletion_folder}: cannot parse rmax/zmax")
            continue

        # Create HDF5 group for this depletion config
        depletion_group = sim_group.require_group(depletion_folder)
        depletion_group.attrs["rmax_kpc"] = rmax
        depletion_group.attrs["zmax_kpc"] = zmax

        # Find the specific depletion_output_XXXXX.txt file
        target_filename = f"depletion_output_{output_number}.txt"
        target_filepath = os.path.join(full_depl_path, target_filename)
        
        if os.path.exists(target_filepath):
            data = load_depletion_file(target_filepath)
            # Create a dataset for this output
            dset_name = f"output_{output_number:05d}"
            if dset_name not in depletion_group:
                depletion_group.create_dataset(dset_name, data=data)
        else:
            print(f"Warning: {target_filename} not found in {depletion_folder}")
        
        # Find the corresponding dust content file
        dust_content_folder = os.path.join(full_depl_path, "DustContent")
        if os.path.exists(dust_content_folder):
            dust_filename = f"dust_content_gal_rcyl{rmax}kpc_zcyl{zmax}kpc_output_{output_number:05d}.txt"
            dust_filepath = os.path.join(dust_content_folder, dust_filename)
            
            if os.path.exists(dust_filepath):
                # Read the single line of dust content data
                with open(dust_filepath, 'r') as f:
                    dust_line = f.readline().strip()
                    dust_values = np.array([float(x) for x in dust_line.split()])

                # Compute the dust-to-metal and dust-to-gas ratios
                dust_values[-2:] = dust_values[-2:] / SioverSil # Correcting the silicate dust masses
                DTM = np.sum(dust_values[-4:]) / dust_values[2] if dust_values[2] != 0 else 0 # Dust-to-metal ratio
                DTG = np.sum(dust_values[-4:]) / dust_values[1] if dust_values[1] != 0 else 0 # Dust-to-gas ratio
                Z   = dust_values[2] / dust_values[1] if dust_values[1] != 0 else 0 # Metallicity
                dust_values = np.append(dust_values, [DTM, DTG, Z])  # Append ratios to the dust values

                # Create a dataset for dust content
                dust_dset_name = f"dust_content_{output_number:05d}"
                if dust_dset_name not in depletion_group:
                    depletion_group.create_dataset(dust_dset_name, data=dust_values)
                    
                    # Add column names as attributes for reference
                    dust_columns = [
                        "time", "gas_mass", "metal_mass", "H_mass", "O_mass", 
                        "Fe_mass", "Mg_mass", "C_mass", "N_mass", "Si_mass", 
                        "S_mass", "D_mass", "Csmall_mass", "Clarge_mass", 
                        "Silsmall_mass", "Sillarge_mass", "DTM", "DTG", "Z"
                    ]
                    depletion_group.attrs["dust_content_columns"] = dust_columns
            else:
                print(f"Warning: {dust_filename} not found in DustContent folder")
        else:
            print(f"Warning: DustContent folder not found in {depletion_folder}")

def read_sim_paths_from_file(file_path):
    """
    Reads a list of simulation paths from a text file.
    Skips empty lines and comments (lines starting with '#').
    """
    with open(file_path, "r") as f:
        lines = f.readlines()
    return [line.strip() for line in lines if line.strip() and not line.startswith("#")]

def main(path_list_file, output_hdf5, output_number):
    """
    Main routine to process all simulation paths and store results in one HDF5 file.
    Only extracts the specified output number.
    """
    sim_paths = read_sim_paths_from_file(path_list_file)

    with h5py.File(output_hdf5, "w") as h5file:
        for sim_path in sim_paths:
            if not os.path.isdir(sim_path):
                print(f"Skipping invalid path: {sim_path}")
                continue
            print(f"Processing simulation: {sim_path}")
            process_simulation(sim_path, h5file, output_number)