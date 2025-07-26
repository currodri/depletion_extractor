import os
import re
import h5py
import numpy as np
from glob import glob

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

def process_simulation(sim_path, h5file):
    """
    Processes one simulation directory:
    - Finds all depletion folders inside it.
    - Parses rmax/zmax from folder names.
    - Extracts all depletion_output_XXXXX.txt files.
    - Stores each output as a dataset inside an HDF5 group.
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

        # Find all depletion_output_XXXXX.txt files
        files = glob(os.path.join(full_depl_path, "depletion_output_*.txt"))
        for filepath in files:
            # Extract output number from the filename
            match = re.search(r"depletion_output_(\d+)\.txt", os.path.basename(filepath))
            if not match:
                continue

            output_number = int(match.group(1))
            data = load_depletion_file(filepath)

            # Create a dataset for this output
            dset_name = f"output_{output_number:05d}"
            if dset_name not in depletion_group:
                depletion_group.create_dataset(dset_name, data=data)

def read_sim_paths_from_file(file_path):
    """
    Reads a list of simulation paths from a text file.
    Skips empty lines and comments (lines starting with '#').
    """
    with open(file_path, "r") as f:
        lines = f.readlines()
    return [line.strip() for line in lines if line.strip() and not line.startswith("#")]

def main(path_list_file, output_hdf5):
    """
    Main routine to process all simulation paths and store results in one HDF5 file.
    """
    sim_paths = read_sim_paths_from_file(path_list_file)

    with h5py.File(output_hdf5, "w") as h5file:
        for sim_path in sim_paths:
            if not os.path.isdir(sim_path):
                print(f"Skipping invalid path: {sim_path}")
                continue
            print(f"Processing simulation: {sim_path}")
            process_simulation(sim_path, h5file)