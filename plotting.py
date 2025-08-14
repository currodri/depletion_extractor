import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import re

# Column names and indices
DEPLETION_NAMES = [
    "rho", "Cdep", "Mgdep", "Fedep", "Sidep", "Odep",
    "varCdep", "varMgdep", "varFedep", "varSidep", "varOdep",
    "C", "Mg", "Fe", "Si", "O", "Cds", "Cdl", "Sds", "Sdl"
]
DEP_IDX = {name: i for i, name in enumerate(DEPLETION_NAMES)}

DUST_CONTENT_NAMES = [
    "time", "gas_mass", "metal_mass", "H_mass", "O_mass",
    "Fe_mass", "Mg_mass", "C_mass", "N_mass", "Si_mass",
    "S_mass", "D_mass", "Csmall_mass", "Clarge_mass",
    "Silsmall_mass", "Sillarge_mass", "DTM", "DTG", "Z"
]
DUST_IDX = {name: i for i, name in enumerate(DUST_CONTENT_NAMES)}

mH_amu = 1.00784 # Hydrogen mass in amu
mO_amu = 15.999 # Oxygen mass in amu

def extract_metallicity_from_name(sim_name):
    """
    Extract initial metallicity from simulation name.
    Expected format: something with 'XXXZsun' where XXX is the metallicity value.
    Returns the metallicity as a float.
    """
    match = re.search(r'(\d*\.?\d+)Zsun', sim_name)
    if match:
        return float(match.group(1))
    else:
        # Default fallback if pattern not found
        return 1.0

def get_simulation_info_from_path(h5_filename):
    """
    Extract simulation type and metallicity information from HDF5 file paths.
    Returns dictionaries mapping sim_name to sim_type and metallicity.
    """
    sim_types = {}
    metallicities = {}
    
    with h5py.File(h5_filename, "r") as f:
        for sim_name in f:
            sim_group = f[sim_name]
            for depletion_config in sim_group:
                group = sim_group[depletion_config]
                
                # Try to get simulation path from HDF5 attributes
                # If path was stored during extraction, use it to determine simulation type
                if 'path' in group.attrs:
                    path = group.attrs['path']
                    # Extract simulation type from path structure: /root/simulation_type/simulation_name
                    path_parts = path.split('/')
                    if len(path_parts) >= 3:
                        # Find the simulation type (e.g., 1d10, 1d11, 1d12) in the path
                        for part in path_parts:
                            if '1d10' in part:
                                sim_types[sim_name] = 'G8'
                                break
                            elif '1d11' in part:
                                sim_types[sim_name] = 'G9'
                                break
                            elif '1d12' in part:
                                sim_types[sim_name] = 'G10'
                                break
                            else:
                                raise ValueError(f"Unknown simulation type in path: {path}")
                    else:
                        raise ValueError(f"Invalid path structure: {path}")
                else:
                    # Fallback: Extract simulation type from simulation name
                    if '1d10' in sim_name:
                        sim_types[sim_name] = 'G8'
                    elif '1d11' in sim_name:
                        sim_types[sim_name] = 'G9'
                    elif '1d12' in sim_name:
                        sim_types[sim_name] = 'G10'
                    else:
                        raise ValueError(f"Unknown simulation type in name: {sim_name}")
                
                # Extract metallicity from name
                metallicity = extract_metallicity_from_name(sim_name)
                metallicities[sim_name] = metallicity
                
                break  # Only need to process one config per simulation
    
    return sim_types, metallicities

def plot_depletion_comparison(h5_filename, output_number):
    """
    Compare element depletions (with error bars) across all simulations
    and depletion regions for a given output number.
    Each subplot shows one element's depletion vs density (nH).
    Uses color coding based on initial metallicity and markers based on simulation type.
    """
    elements = ["Cdep", "Mgdep", "Fedep", "Sidep", "Odep"]
    var_elements = ["varCdep", "varMgdep", "varFedep", "varSidep", "varOdep"]
    n_elements = len(elements)

    # Get simulation information
    sim_types, metallicities = get_simulation_info_from_path(h5_filename)
    print(f"Found {len(sim_types)} simulations with types and metallicities.")
    print(f"Metallicities: {metallicities}")
    print(f"Simulations: {list(sim_types.keys())}")
    
    # Define markers and colormaps for different simulation types
    type_config = {
        'G8': {'marker': 'o', 'cmap': cm.Greens},
        'G9': {'marker': 's', 'cmap': cm.Blues}, 
        'G10': {'marker': '^', 'cmap': cm.Purples},
        'Unknown': {'marker': 'D', 'cmap': cm.Greys}
    }
    
    # Get unique simulation types and assign markers/colormaps
    unique_types = list(set(sim_types.values()))
    available_markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
    available_cmaps = [cm.Greens, cm.Blues, cm.Purples, cm.Greys, cm.Oranges, 
                       cm.Reds, cm.YlOrBr, cm.YlOrRd, cm.BuPu, cm.GnBu]
    
    for i, sim_type in enumerate(unique_types):
        if sim_type not in type_config:
            type_config[sim_type] = {
                'marker': available_markers[i % len(available_markers)],
                'cmap': available_cmaps[i % len(available_cmaps)]
            }
    
    # Set up normalization for metallicity values
    metallicity_values = list(metallicities.values())
    if len(set(metallicity_values)) > 1:
        norm = plt.Normalize(vmin=min(metallicity_values), vmax=max(metallicity_values))
    else:
        # If all metallicities are the same, use a single color
        norm = plt.Normalize(vmin=0, vmax=1)

    fig, axes = plt.subplots(n_elements, 1, figsize=(5, 3* n_elements), sharex=True)

    with h5py.File(h5_filename, "r") as f:
        for sim_name in f:
            sim_group = f[sim_name]
            for depletion_config in sim_group:
                group = sim_group[depletion_config]
                output_key = f"output_{int(output_number):05d}"
                if output_key not in group:
                    continue

                # Get marker, colormap, and color for this simulation
                sim_type = sim_types[sim_name]
                marker = type_config[sim_type]['marker']
                cmap = type_config[sim_type]['cmap']
                metallicity = metallicities[sim_name]
                color = cmap(norm(metallicity))
                
                data = group[output_key][()]
                rho = data[:, DEP_IDX["rho"]]

                for i, (elem, var) in enumerate(zip(elements, var_elements)):
                    dep = data[:, DEP_IDX[elem]]
                    err = np.sqrt(data[:, DEP_IDX[var]])
                    axes[i].errorbar(
                        rho, dep, yerr=err, fmt=f"{marker}-", markersize=10,
                        alpha=0.7, capsize=4, color=color
                    )

    # Create custom legend showing simulation types only
    legend_elements = []
    unique_sim_types = list(set(sim_types.values()))
    for sim_type in unique_sim_types:
        marker = type_config[sim_type]['marker']
        cmap = type_config[sim_type]['cmap']
        # Use middle color of the colormap as representative
        representative_color = cmap(0.5)
        legend_elements.append(plt.Line2D([0], [0], marker=marker, color=representative_color,
                                        linestyle='-', markersize=10, label=f'{sim_type}'))

    # Format plots
    for i, elem in enumerate(elements):
        ax = axes[i]
        ax.set_ylabel(f"{elem} (fraction)")
        ax.set_xscale("log")
        ax.set_ylim(0, 1)
        ax.grid(True, which="both", ls="--", lw=0.3)
        # Add also minor ticks
        ax.minorticks_on()
        ax.tick_params(which='both', direction='in', top=True, right=True)
        if i == 0:  # Only add legend to the first subplot
            ax.legend(handles=legend_elements, fontsize=10, loc="lower left", frameon=False)

    axes[-1].set_xlabel("Density (nH)")
    fig.suptitle(f"Depletion Comparison at Output {int(output_number):05d}", fontsize=15)
    
    plt.tight_layout()
    plt.savefig(f"depletion_comparison_output_{int(output_number):05d}.png", format='png', dpi=300)

def plot_dtm_dtg_vs_metallicity(h5_filename, output_number):
    """
    Plot dust-to-metal ratio (DTM) vs metallicity for each simulation
    and depletion region at a given output number.
    Uses color coding based on initial metallicity and markers based on simulation type.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4))

    OH_Asplund2009  = 4.9e-4 # Solar oxygen abundance (Asplund et al. 2009)

    # Get simulation information
    sim_types, metallicities = get_simulation_info_from_path(h5_filename)
    
    # Define markers and colormaps for different simulation types
    type_config = {
        'G8': {'marker': 'o', 'cmap': cm.Greens},
        'G9': {'marker': 's', 'cmap': cm.Blues}, 
        'G10': {'marker': '^', 'cmap': cm.Purples},
        'Unknown': {'marker': 'D', 'cmap': cm.Greys}
    }
    
    # Get unique simulation types and assign markers/colormaps
    unique_types = list(set(sim_types.values()))
    available_markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
    available_cmaps = [cm.Greens, cm.Blues, cm.Purples, cm.Greys, cm.Oranges, 
                       cm.Reds, cm.YlOrBr, cm.YlOrRd, cm.BuPu, cm.GnBu]
    
    for i, sim_type in enumerate(unique_types):
        if sim_type not in type_config:
            type_config[sim_type] = {
                'marker': available_markers[i % len(available_markers)],
                'cmap': available_cmaps[i % len(available_cmaps)]
            }
    
    # Set up normalization for metallicity values
    metallicity_values = list(metallicities.values())
    if len(set(metallicity_values)) > 1:
        norm = plt.Normalize(vmin=min(metallicity_values), vmax=max(metallicity_values))
    else:
        # If all metallicities are the same, use a single color
        norm = plt.Normalize(vmin=0, vmax=1)

    with h5py.File(h5_filename, "r") as f:
        for sim_name in f:
            sim_group = f[sim_name]
            for depletion_config in sim_group:
                group = sim_group[depletion_config]
                output_key = f"dust_content_{int(output_number):05d}"
                if output_key not in group:
                    continue

                # Get marker, colormap, and color for this simulation
                sim_type = sim_types[sim_name]
                marker = type_config[sim_type]['marker']
                cmap = type_config[sim_type]['cmap']
                metallicity = metallicities[sim_name]
                color = cmap(norm(metallicity))
                
                data = group[output_key][()]

                # Convert Oxygen and Hydrogen masses to abundances
                OH = data[DUST_IDX["O_mass"]] * mH_amu / (data[DUST_IDX["H_mass"]] * mO_amu)
                metallicity_measured = OH / OH_Asplund2009  # Relative to solar abundance
                dtm = data[DUST_IDX["DTM"]]
                dtg = data[DUST_IDX["DTG"]]

                ax1.scatter(metallicity_measured, dtm, marker=marker, color=color, s=100)
                ax2.scatter(metallicity_measured, dtg, marker=marker, color=color, s=100)

    # Create custom legend showing simulation types only
    legend_elements = []
    unique_sim_types = list(set(sim_types.values()))
    for sim_type in unique_sim_types:
        marker = type_config[sim_type]['marker']
        cmap = type_config[sim_type]['cmap']
        # Use middle color of the colormap as representative
        representative_color = cmap(0.5)
        legend_elements.append(plt.Line2D([0], [0], marker=marker, color=representative_color,
                                        linestyle='None', markersize=12, label=f'{sim_type}'))

    ax1.set_xlabel(r"$(O/H)/(O/H)_\odot$")
    ax1.set_ylabel("Dust-to-Metal Ratio (DTM)")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid(True, which="both", ls="--", lw=0.3)
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)
    ax1.legend(handles=legend_elements, fontsize=10, loc="upper left", frameon=False)
    
    ax2.set_xlabel(r"$(O/H)/(O/H)_\odot$")
    ax2.set_ylabel("Dust-to-Gas Ratio (DTG)")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.grid(True, which="both", ls="--", lw=0.3)
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)
    # No legend on second subplot to avoid duplication
    
    fig.suptitle(f"DTM and DTG vs Metallicity at Output {int(output_number):05d}", fontsize=15)
    
    plt.tight_layout()
    plt.savefig(f"dtm_dtg_vs_metallicity_output_{int(output_number):05d}.png", format='png', dpi=300)