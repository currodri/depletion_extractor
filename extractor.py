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
    Extracts the last two components of the simulation path.
    E.g., '/data/simtype/sim1' -> 'simtype/sim1'
    """
    abs_path = os.path.abspath(path)
    parts = abs_path.strip(os.sep).split(os.sep)
    if len(parts) >= 2:
        return os.path.join(parts[-2], parts[-1])
    else:
        return parts[-1] if parts else ""

def load_depletion_file(filepath):
    """
    Loads a depletion file into a 2D NumPy array.
    Assumes the file is plain text with space-separated values.
    """
    return np.loadtxt(filepath)

def process_simulation(sim_path, h5file, output_number, target_rmax_depl, target_zmax_depl, target_rmax_dust, target_zmax_dust):
    """
    Processes one simulation directory for specific rmax/zmax configurations:
    - Finds the specific depletion folder matching rmax_depl/zmax_depl.
    - Extracts only the specified depletion_output_XXXXX.txt file.
    - Finds dust content file matching rmax_dust/zmax_dust.
    - Stores the output as a dataset inside an HDF5 group.
    """
    sim_name = get_sim_name(sim_path)
    sim_group = h5file.require_group(sim_name)

    # Construct the expected depletion folder name
    depletion_folder = f"DEPLETION_rmax{target_rmax_depl}kpc_zmax{target_zmax_depl}kpc"
    full_depl_path = os.path.join(sim_path, depletion_folder)
    
    if not os.path.exists(full_depl_path):
        print(f"Warning: {depletion_folder} not found in {sim_path}")
        return

    # Create HDF5 group for this depletion config
    depletion_group = sim_group.require_group(depletion_folder)
    depletion_group.attrs["rmax_depl_kpc"] = target_rmax_depl
    depletion_group.attrs["zmax_depl_kpc"] = target_zmax_depl
    depletion_group.attrs["rmax_dust_kpc"] = target_rmax_dust
    depletion_group.attrs["zmax_dust_kpc"] = target_zmax_dust

    # Find the specific depletion_output_XXXXX.txt file
    target_filename = f"depletion_output_{output_number:05d}.txt"
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
    dust_content_folder = os.path.join(sim_path, "DustContent")
    if os.path.exists(dust_content_folder):
        dust_filename = f"dust_content_gal_rcyl{target_rmax_dust}kpc_zcyl{target_zmax_dust}kpc_output_{output_number:05d}.txt"
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
            print(f"DTM: {DTM}, DTG: {DTG} for {dust_filename}")
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
        print(f"Warning: DustContent folder not found in {sim_path}")

def read_sim_paths_from_file(file_path):
    """
    Reads a list of simulation paths with rmax and zmax from a text file.
    Expected format: path rmax_depl zmax_depl rmax_dust zmax_dust (space-separated)
    Skips empty lines and comments (lines starting with '#').
    Returns: List of tuples (path, rmax_depl, zmax_depl, rmax_dust, zmax_dust)
    """
    sim_configs = []
    with open(file_path, "r") as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        
        parts = line.split()
        if len(parts) != 5:
            print(f"Warning: Invalid format in line '{line}'. Expected: path rmax_depl zmax_depl rmax_dust zmax_dust")
            continue
        
        try:
            path = parts[0]
            rmax_depl = float(parts[1])
            zmax_depl = float(parts[2])
            rmax_dust = float(parts[3])
            zmax_dust = float(parts[4])
            sim_configs.append((path, rmax_depl, zmax_depl, rmax_dust, zmax_dust))
        except ValueError:
            print(f"Warning: Could not parse rmax/zmax values from line '{line}'")
            continue
    
    return sim_configs

def main(path_list_file, output_hdf5, output_number):
    """
    Main routine to process all simulation paths with specific rmax/zmax and store results in one HDF5 file.
    Only extracts the specified output number for the specified depletion configurations.
    """
    sim_configs = read_sim_paths_from_file(path_list_file)

    with h5py.File(output_hdf5, "w") as h5file:
        for sim_path, rmax_depl, zmax_depl, rmax_dust, zmax_dust in sim_configs:
            if not os.path.isdir(sim_path):
                print(f"Skipping invalid path: {sim_path}")
                continue
            print(f"Processing simulation: {sim_path} (depl: rmax={rmax_depl}, zmax={zmax_depl}; dust: rmax={rmax_dust}, zmax={zmax_dust})")
            process_simulation(sim_path, h5file, output_number, rmax_depl, zmax_depl, rmax_dust, zmax_dust)