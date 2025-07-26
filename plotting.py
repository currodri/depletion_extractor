import h5py
import numpy as np
import matplotlib.pyplot as plt

# Column names and indices
COLUMN_NAMES = [
    "rho", "Cdep", "Mgdep", "Fedep", "Sidep", "Odep",
    "varCdep", "varMgdep", "varFedep", "varSidep", "varOdep",
    "C", "Mg", "Fe", "Si", "O", "Cds", "Cdl", "Sds", "Sdl"
]
COL_IDX = {name: i for i, name in enumerate(COLUMN_NAMES)}

def plot_depletion_comparison(h5_filename, output_number):
    """
    Compare element depletions (with error bars) across all simulations
    and depletion regions for a given output number.
    Each subplot shows one element's depletion vs density (nH).
    """
    elements = ["Cdep", "Mgdep", "Fedep", "Sidep", "Odep"]
    var_elements = ["varCdep", "varMgdep", "varFedep", "varSidep", "varOdep"]
    n_elements = len(elements)

    fig, axes = plt.subplots(n_elements, 1, figsize=(9, 4 * n_elements), sharex=True)

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
                rho = data[:, COL_IDX["rho"]]

                for i, (elem, var) in enumerate(zip(elements, var_elements)):
                    dep = data[:, COL_IDX[elem]]
                    err = np.sqrt(data[:, COL_IDX[var]])
                    axes[i].errorbar(
                        rho, dep, yerr=err, fmt="o-", markersize=4,
                        label=label, alpha=0.7, capsize=2
                    )

    # Format plots
    for i, elem in enumerate(elements):
        ax = axes[i]
        ax.set_ylabel(f"{elem} (fraction)")
        ax.set_xscale("log")
        ax.set_ylim(0, 1)
        ax.grid(True, which="both", ls="--", lw=0.3)
        ax.legend(fontsize=8, loc="lower left")

    axes[-1].set_xlabel("Density (nH)")
    fig.suptitle(f"Depletion Comparison at Output {int(output_number):05d}", fontsize=15)
    plt.tight_layout()
    plt.savefig(f"depletion_comparison_output_{int(output_number):05d}.png")