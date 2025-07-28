import h5py
import numpy as np
import matplotlib.pyplot as plt

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

def plot_depletion_comparison(h5_filename, output_number):
    """
    Compare element depletions (with error bars) across all simulations
    and depletion regions for a given output number.
    Each subplot shows one element's depletion vs density (nH).
    """
    elements = ["Cdep", "Mgdep", "Fedep", "Sidep", "Odep"]
    var_elements = ["varCdep", "varMgdep", "varFedep", "varSidep", "varOdep"]
    n_elements = len(elements)

    fig, axes = plt.subplots(n_elements, 1, figsize=(5, 3* n_elements), sharex=True)

    with h5py.File(h5_filename, "r") as f:
        for sim_name in f:
            sim_group = f[sim_name]
            for depletion_config in sim_group:
                group = sim_group[depletion_config]
                output_key = f"output_{int(output_number):05d}"
                if output_key not in group:
                    continue

                label = f"{sim_name} / {depletion_config}"
                data = group[output_key][()]
                rho = data[:, DEP_IDX["rho"]]

                for i, (elem, var) in enumerate(zip(elements, var_elements)):
                    dep = data[:, DEP_IDX[elem]]
                    err = np.sqrt(data[:, DEP_IDX[var]])
                    axes[i].errorbar(
                        rho, dep, yerr=err, fmt="o-", markersize=10,
                        label=label, alpha=0.7, capsize=4
                    )

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
        ax.legend(fontsize=8, loc="lower left",frameon=False)

    axes[-1].set_xlabel("Density (nH)")
    fig.suptitle(f"Depletion Comparison at Output {int(output_number):05d}", fontsize=15)
    plt.tight_layout()
    plt.savefig(f"depletion_comparison_output_{int(output_number):05d}.png",format='png', dpi=300)

def plot_dtm_dtg_vs_metallicity(h5_filename, output_number):
    """
    Plot dust-to-metal ratio (DTM) vs metallicity for each simulation
    and depletion region at a given output number.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    OH_Asplund2009  = 8.69  # Solar oxygen abundance (Asplund et al. 2009)

    with h5py.File(h5_filename, "r") as f:
        for sim_name in f:
            sim_group = f[sim_name]
            for depletion_config in sim_group:
                group = sim_group[depletion_config]
                output_key = f"output_{int(output_number):05d}"
                if output_key not in group:
                    continue

                label = f"{sim_name} / {depletion_config}"
                data = group[output_key][()]
                # Convert Oxygen and Hydrogen masses to abundances
                OH = data[:, DUST_IDX["O_mass"]] * mH_amu / (data[:, DUST_IDX["H_mass"]] * mO_amu)
                metallicity = OH / OH_Asplund2009  # Relative to solar abundance
                dtm = data[:, DUST_IDX["DTM"]]
                dtg = data[:, DUST_IDX["DTG"]]

                ax1.plot(metallicity, dtm, 'o-', label=label)
                ax2.plot(metallicity, dtg, 'o-', label=label)

    ax1.set_xlabel(r"$(O/H)/(O/H)_\odot$")
    ax1.set_ylabel("Dust-to-Metal Ratio (DTM)")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid(True, which="both", ls="--", lw=0.3)
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)
    ax1.legend(fontsize=8, loc="upper left", frameon=False)
    ax2.set_xlabel(r"$(O/H)/(O/H)_\odot$")
    ax2.set_ylabel("Dust-to-Gas Ratio (DTG)")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.grid(True, which="both", ls="--", lw=0.3)
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)
    ax2.legend(fontsize=8, loc="upper left", frameon=False)
    fig.suptitle(f"DTM and DTG vs Metallicity at Output {int(output_number):05d}", fontsize=15)
    plt.tight_layout()
    plt.savefig(f"dtm_dtg_vs_metallicity_output_{int(output_number):05d}.png", format='png', dpi=300)