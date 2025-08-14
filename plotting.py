import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
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
    Extract simulation type and metallicity information from HDF5 file attributes.
    Returns dictionaries mapping sim_name to sim_type and metallicity.
    """
    sim_types = {}
    metallicities = {}
    
    with h5py.File(h5_filename, "r") as f:
        for sim_name in f:
            sim_group = f[sim_name]
            
            # Check if simulation metadata is stored as attributes
            if 'simulation_type' in sim_group.attrs and 'initial_metallicity' in sim_group.attrs:
                # Use pre-stored attributes (preferred method)
                sim_type = sim_group.attrs['simulation_type']
                if sim_type not in ['G8', 'G9', 'G10']:
                    raise ValueError(f"Invalid simulation type '{sim_type}' for simulation '{sim_name}'. Must be one of G8, G9, or G10.")
                sim_types[sim_name] = sim_type
                metallicities[sim_name] = sim_group.attrs['initial_metallicity']
            else:
                # Fallback: Extract from name and path (for backward compatibility)
                sim_type = None
                for depletion_config in sim_group:
                    group = sim_group[depletion_config]
                    
                    # Try to get simulation path from HDF5 attributes
                    if 'path' in group.attrs:
                        path = group.attrs['path']
                        # Extract simulation type from path structure
                        path_parts = path.split('/')
                        for part in path_parts:
                            if part in ['G8', 'G9', 'G10']:
                                sim_type = part
                                break
                    
                    # If not found in path, try simulation name
                    if sim_type is None:
                        if 'G8' in sim_name or '1d10' in sim_name:
                            sim_type = 'G8'
                        elif 'G9' in sim_name or '1d11' in sim_name:
                            sim_type = 'G9'
                        elif 'G10' in sim_name or '1d12' in sim_name:
                            sim_type = 'G10'
                    
                    # Raise error if simulation type cannot be determined
                    if sim_type is None:
                        raise ValueError(f"Cannot determine simulation type for simulation '{sim_name}'. Must be one of G8, G9, or G10.")
                    
                    sim_types[sim_name] = sim_type
                    
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
    
    # Define colormaps for different simulation types (markers determined by gas fraction)
    type_config = {
        'G8': {'cmap': cm.Greens},
        'G9': {'cmap': cm.Blues}, 
        'G10': {'cmap': cm.Purples}
    }
    
    # Get unique simulation types and verify they are supported
    unique_types = list(set(sim_types.values()))
    for sim_type in unique_types:
        if sim_type not in type_config:
            raise ValueError(f"Unsupported simulation type '{sim_type}'. Only G8, G9, and G10 are supported.")
    
    # Set up normalization for metallicity values per simulation type (in log scale)
    type_norms = {}
    for sim_type in unique_types:
        # Get metallicities for this simulation type only
        type_metallicities = [metallicities[sim_name] for sim_name, s_type in sim_types.items() if s_type == sim_type]
        if len(set(type_metallicities)) > 1:
            type_norms[sim_type] = LogNorm(vmin=min(type_metallicities), vmax=max(type_metallicities))
        else:
            # If all metallicities are the same for this type, use a single color
            type_norms[sim_type] = plt.Normalize(vmin=0, vmax=1)

    fig, axes = plt.subplots(n_elements, 1, figsize=(5, 3* n_elements), sharex=True)

    with h5py.File(h5_filename, "r") as f:
        for sim_name in f:
            sim_group = f[sim_name]
            for depletion_config in sim_group:
                group = sim_group[depletion_config]
                output_key = f"output_{int(output_number):05d}"
                if output_key not in group:
                    continue

                # Get marker based on gas fraction type and colormap/color for this simulation
                sim_type = sim_types[sim_name]
                # Determine marker based on gas fraction: 'fg' in name = low gas fraction (squares), else regular (circles)
                marker = 's' if 'fg' in sim_name.lower() else 'o'
                cmap = type_config[sim_type]['cmap']
                metallicity = metallicities[sim_name]
                norm = type_norms[sim_type]  # Use the normalization for this simulation type
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

    # Create custom legend showing simulation types and gas fraction types
    legend_elements = []
    unique_sim_types = list(set(sim_types.values()))
    
    # Add simulation type legend entries (using representative colors)
    for sim_type in unique_sim_types:
        cmap = type_config[sim_type]['cmap']
        # Use middle color of the colormap as representative
        representative_color = cmap(0.5)
        legend_elements.append(plt.Line2D([0], [0], marker='o', color=representative_color,
                                        linestyle='-', markersize=10, label=f'{sim_type}'))
    
    # Add gas fraction type legend entries
    legend_elements.append(plt.Line2D([0], [0], marker='o', color='gray',
                                    linestyle='-', markersize=10, label='Regular gas fraction'))
    legend_elements.append(plt.Line2D([0], [0], marker='s', color='gray',
                                    linestyle='-', markersize=10, label='Low gas fraction (fg)'))

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
    
    # Define colormaps for different simulation types (markers determined by gas fraction)
    type_config = {
        'G8': {'cmap': cm.Greens},
        'G9': {'cmap': cm.Blues}, 
        'G10': {'cmap': cm.Purples}
    }
    
    # Get unique simulation types and verify they are supported
    unique_types = list(set(sim_types.values()))
    for sim_type in unique_types:
        if sim_type not in type_config:
            raise ValueError(f"Unsupported simulation type '{sim_type}'. Only G8, G9, and G10 are supported.")
    
    # Set up normalization for metallicity values per simulation type (in log scale)
    type_norms = {}
    for sim_type in unique_types:
        # Get metallicities for this simulation type only
        type_metallicities = [metallicities[sim_name] for sim_name, s_type in sim_types.items() if s_type == sim_type]
        if len(set(type_metallicities)) > 1:
            # Use log normalization
            type_norms[sim_type] = LogNorm(vmin=min(type_metallicities), vmax=max(type_metallicities))
        else:
            # If all metallicities are the same for this type, use a single color
            type_norms[sim_type] = plt.Normalize(vmin=0, vmax=1)

    with h5py.File(h5_filename, "r") as f:
        for sim_name in f:
            sim_group = f[sim_name]
            for depletion_config in sim_group:
                group = sim_group[depletion_config]
                output_key = f"dust_content_{int(output_number):05d}"
                if output_key not in group:
                    continue

                # Get marker based on gas fraction type and colormap/color for this simulation
                sim_type = sim_types[sim_name]
                # Determine marker based on gas fraction: 'fg' in name = low gas fraction (squares), else regular (circles)
                marker = 's' if 'fg' in sim_name.lower() else 'o'
                cmap = type_config[sim_type]['cmap']
                metallicity = metallicities[sim_name]
                norm = type_norms[sim_type]  # Use the normalization for this simulation type
                color = cmap(norm(metallicity))
                
                data = group[output_key][()]

                # Convert Oxygen and Hydrogen masses to abundances
                OH = data[DUST_IDX["O_mass"]] * mH_amu / (data[DUST_IDX["H_mass"]] * mO_amu)
                metallicity_measured = OH / OH_Asplund2009  # Relative to solar abundance
                dtm = data[DUST_IDX["DTM"]]
                dtg = data[DUST_IDX["DTG"]]

                ax1.scatter(metallicity_measured, dtm, marker=marker, color=color, s=100)
                ax2.scatter(metallicity_measured, dtg, marker=marker, color=color, s=100)

    # Create custom legend showing simulation types and gas fraction types
    legend_elements = []
    unique_sim_types = list(set(sim_types.values()))
    
    # Add simulation type legend entries (using representative colors)
    for sim_type in unique_sim_types:
        cmap = type_config[sim_type]['cmap']
        # Use middle color of the colormap as representative
        representative_color = cmap(0.5)
        legend_elements.append(plt.Line2D([0], [0], marker='o', color=representative_color,
                                        linestyle='None', markersize=12, label=f'{sim_type}'))
    
    # Add gas fraction type legend entries
    legend_elements.append(plt.Line2D([0], [0], marker='o', color='gray',
                                    linestyle='None', markersize=12, label='Regular gas fraction'))
    legend_elements.append(plt.Line2D([0], [0], marker='s', color='gray',
                                    linestyle='None', markersize=12, label='Low gas fraction (fg)'))

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