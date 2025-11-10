<div align="center" id="top"> 
  <img src="./.github/app.gif" alt="Depletion_extractor" />

  &#xa0;

  <!-- <a href="https://depletion_extractor.netlify.app">Demo</a> -->
</div>

<h1 align="center">Depletion Extractor</h1>

<p align="center">
  <img alt="Github top language" src="https://img.shields.io/github/languages/top/currodri/depletion_extractor?color=56BEB8">

  <img alt="Github language count" src="https://img.shields.io/github/languages/count/currodri/depletion_extractor?color=56BEB8">

  <img alt="Repository size" src="https://img.shields.io/github/repo-size/currodri/depletion_extractor?color=56BEB8">

  <img alt="License" src="https://img.shields.io/github/license/currodri/depletion_extractor?color=56BEB8">

  <!-- <img alt="Github issues" src="https://img.shields.io/github/issues/currodri/depletion_extractor?color=56BEB8" /> -->

  <!-- <img alt="Github forks" src="https://img.shields.io/github/forks/currodri/depletion_extractor?color=56BEB8" /> -->

  <!-- <img alt="Github stars" src="https://img.shields.io/github/stars/currodri/depletion_extractor?color=56BEB8" /> -->
</p>

<!-- Status -->

<!-- <h4 align="center"> 
	🚧  Depletion_extractor 🚀 Under construction...  🚧
</h4> 

<hr> -->

<p align="center">
  <a href="#dart-about">About</a> &#xa0; | &#xa0; 
  <a href="#sparkles-features">Features</a> &#xa0; | &#xa0;
  <a href="#rocket-technologies">Technologies</a> &#xa0; | &#xa0;
  <a href="#white_check_mark-requirements">Requirements</a> &#xa0; | &#xa0;
  <a href="#checkered_flag-starting">Starting</a> &#xa0; | &#xa0;
  <a href="#memo-license">License</a> &#xa0; | &#xa0;
  <a href="https://github.com/currodri" target="_blank">Author</a>
</p>

<br>

## :dart: About ##

The tools here included allow the extraction of average dust and metal analysis of the suite of G8, G9, and G10 isolated galaxy simulations presented in Dubois et al. (2024) ([paper link](https://www.aanda.org/articles/aa/full_html/2024/07/aa49784-24/aa49784-24.html)). These simulations have been run with the adaptive mesh-refinement code for hydrodynamics and gravity RAMSES, modified to include on-the-fly modelling of dust evolution developed by Yohan Dubois (IAP) and Curro Rodriguez Montero (KICP, UChicago).

If you use this dataset, please cite the original work of Dubois et al. (2024) ([paper link](https://www.aanda.org/articles/aa/full_html/2024/07/aa49784-24/aa49784-24.html)). For any questions about this code, the dataset or the simulations, please feel free to contact Curro Rodriguez Montero ([currodri@uchicago.edu](mailto:currodri@uchicago.edu)) or Yohan Dubois ([dubois@iap.fr](mailto:dubois@iap.fr)).

## :sparkles: Features ##

:heavy_check_mark: Extract element depletion data from RAMSES simulation outputs;\
:heavy_check_mark: Process dust content and mass data for different galactic regions;\
:heavy_check_mark: Calculate dust-to-metal and dust-to-gas ratios;\
:heavy_check_mark: Generate comparison plots across simulations;\
:heavy_check_mark: Efficient HDF5 storage format for large datasets;

## :rocket: Technologies ##

The following tools were used in this project:

- [Python](https://python.org/) - Main programming language
- [NumPy](https://numpy.org/) - Numerical computations and array handling
- [h5py](https://www.h5py.org/) - HDF5 file format handling
- [Matplotlib](https://matplotlib.org/) - Data visualization and plotting
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) - Efficient binary data storage format

## :white_check_mark: Requirements ##

Before starting :checkered_flag:, you need to have [Python 3.6+](https://python.org/) installed with the following packages:

```bash
pip install numpy h5py matplotlib
```

## :package: What does this package do? ##

This package extracts and processes **element depletion data** from cosmological galaxy simulations. Element depletion refers to the process where metals (elements heavier than hydrogen and helium) are removed from the gas phase and incorporated into dust grains in the interstellar medium.

### Sample Dataset Included

The package includes a sample HDF5 dataset (`sample_data.h5`) containing processed simulation data that you can use immediately to test all functionality and reproduce the examples in this README. The dataset is small (95 KB) and contains representative data from multiple G8, G9, and G10 simulations with different metallicities and gas fraction configurations.

### Key Functionality:

1. **Data Extraction**: Reads depletion output files from RAMSES simulations for specific spatial regions defined by radial (`rmax`) and vertical (`zmax`) extents
2. **Dust Content Analysis**: Processes dust mass data for different dust species (carbon small/large grains, silicate small/large grains)
3. **Ratio Calculations**: Computes important astrophysical ratios including:
   - Dust-to-metal ratio (DTM)
   - Dust-to-gas ratio (DTG) 
   - Metallicity (Z)
4. **Data Organization**: Stores all data in efficient HDF5 format for fast access and analysis
5. **Visualization**: Generates comparison plots showing depletion trends across simulations

## :file_folder: HDF5 File Structure ##

The output HDF5 file is organized hierarchically as follows:

```
output_file.h5
├── simulation_type_simulation_name/          # Simulation names now use underscore separator
│   ├── attributes:                           # Simulation-level metadata
│   │   ├── simulation_type                   # G8, G9, or G10
│   │   ├── initial_metallicity               # Initial metallicity value
│   │   └── original_path                     # Full path to simulation directory
│   └── DEPLETION_rmax{X}kpc_zmax{Y}kpc/
│       ├── output_{NNNNN}                    # 2D array: depletion data vs density
│       ├── dust_content_{NNNNN}              # 1D array: integrated dust/gas properties
│       └── attributes:
│           ├── rmax_depl_kpc                 # Radial extent for depletion analysis
│           ├── zmax_depl_kpc                 # Vertical extent for depletion analysis
│           ├── rmax_dust_kpc                 # Radial extent for dust content analysis
│           ├── zmax_dust_kpc                 # Vertical extent for dust content analysis
│           ├── path                          # Original simulation path
│           └── dust_content_columns          # Column names for dust_content data
├── simulation_type_simulation_name/
│   └── DEPLETION_rmax{X}kpc_zmax{Y}kpc/
│       └── ...
└── ...
```

### Depletion Data (`output_{NNNNN}`)
2D array with columns:
- `rho`: Gas density (nH in cm⁻³)
- `Cdep`, `Mgdep`, `Fedep`, `Sidep`, `Odep`: Depletion fractions for C, Mg, Fe, Si, O
- `varCdep`, `varMgdep`, `varFedep`, `varSidep`, `varOdep`: Variance in depletion measurements
- `C`, `Mg`, `Fe`, `Si`, `O`: Absolute abundances
- `Cds`, `Cdl`, `Sds`, `Sdl`: Small/large grain dust species masses

### Dust Content Data (`dust_content_{NNNNN}`)
1D array with integrated quantities:
- `time`: Simulation time
- `gas_mass`: Total gas mass in region
- `metal_mass`: Total metal mass in region  
- `H_mass` through `D_mass`: Individual element masses
- `Csmall_mass`, `Clarge_mass`: Carbon dust grain masses
- `Silsmall_mass`, `Sillarge_mass`: Silicate dust grain masses (corrected for Si/silicate ratio)
- `DTM`: Dust-to-metal ratio
- `DTG`: Dust-to-gas ratio
- `Z`: Gas-phase metallicity

## :checkered_flag: Usage ##

### Quick Start with Sample Data

The easiest way to get started is using the included sample dataset:

```bash
# Clone this project
$ git clone https://github.com/currodri/depletion_extractor

# Access the directory
$ cd depletion_extractor

# Test the plotting functions with the sample data
$ python -c "
from plotting import plot_depletion_comparison, plot_dtm_dtg_vs_metallicity
plot_depletion_comparison('sample_data.h5', output_number=81)
plot_dtm_dtg_vs_metallicity('sample_data.h5', output_number=81)
"

# Or run the examples interactively
$ python
>>> from plotting import plot_depletion_comparison, plot_dtm_dtg_vs_metallicity
>>> plot_depletion_comparison('sample_data.h5', 81)
>>> plot_dtm_dtg_vs_metallicity('sample_data.h5', 81)
```

This will generate comparison plots using the sample data, demonstrating the visualization capabilities without needing to process your own simulation data.

### Processing Your Own Data

To extract data from your own RAMSES simulations:

### Input File Format

Create a text file listing simulation paths with their corresponding rmax/zmax values for both depletion and dust content analysis:

```
# path_list.txt
/path/to/simulation1 4.0 1.0 4.0 1.0
/path/to/simulation2 8.0 0.2 8.0 0.2  
/path/to/simulation3 12.0 0.5 12.0 0.5
```

Each line contains: `simulation_path rmax_depl zmax_depl rmax_dust zmax_dust`

Note: The tool now supports different rmax/zmax values for depletion and dust content extraction, allowing for more flexible analysis of different galactic regions.

### Running the Extraction

```bash
# Clone this project
$ git clone https://github.com/currodri/depletion_extractor

# Access the directory
$ cd depletion_extractor

# Run the extraction and plotting
$ python main.py path_list.txt -f output_data.h5 -o 81

# Arguments:
# path_list.txt     : File containing simulation paths and parameters
# -f output_data.h5 : Output HDF5 filename (optional, default: sample_data.h5)
# -o 81            : Output number to extract (optional, default: 81)
```

### Output specification (single, list, or range)

The `-o/--outputs` option now accepts either a single output number, a comma-separated list, or an inclusive range. Examples:

```bash
# Single output (backwards compatible)
$ python main.py path_list.txt -f output_data.h5 -o 81

# Comma-separated list of specific outputs
$ python main.py path_list.txt -f output_data.h5 -o 10,20,30

# Inclusive range of outputs (extracts every output between 10 and 20)
$ python main.py path_list.txt -f output_data.h5 -o 10-20

# Combined list+range (mix works too)
$ python main.py path_list.txt -f output_data.h5 -o 5,10-12,20
```

Notes:
- When multiple outputs are requested, the extractor will write datasets for each output into the same HDF5 file under each simulation's depletion group (e.g. `output_00010`, `output_00011`, ...).
- Automatic plotting is only performed when a single output is requested (to avoid clutter). If you supply multiple outputs, the script will skip plotting and only produce the HDF5 file; you can still plot any extracted output later by calling the plotting functions or re-running `main.py` with a single `-o` value.

### Expected Directory Structure

Each simulation directory should contain:
```
simulation_directory/
├── DEPLETION_rmax{X}kpc_zmax{Y}kpc/
│   ├── depletion_output_{NNNNN}.txt    # Depletion vs density data
│   └── DustContent/
│       └── dust_content_gal_rcyl{X}kpc_zcyl{Y}kpc_output_{NNNNN}.txt
└── DustContent/                        # Alternative location for dust content files
    └── dust_content_gal_rcyl{X}kpc_zcyl{Y}kpc_output_{NNNNN}.txt
```

The tool will automatically generate:
- HDF5 file with extracted data
- Comparison plots (PNG format) showing depletion trends

## :chart_with_upwards_trend: Advanced Plotting and Visualization ##

The package includes sophisticated plotting functions that use color coding and marker shapes to distinguish between simulation types, metallicities, and gas fraction models:

### Built-in Plotting Functions

The package provides two main plotting functions:

1. **`plot_depletion_comparison()`**: Shows element depletion vs density for all simulations
2. **`plot_dtm_dtg_vs_metallicity()`**: Displays dust-to-metal and dust-to-gas ratios vs measured metallicity

### Visual Encoding System

- **Colors**: Different simulation types (G8=Greens, G9=Blues, G10=Purples)
- **Color intensity**: Initial metallicity within each simulation type (log scale normalization)
- **Marker shapes**: 
  - Circles (○): Regular gas fraction simulations
  - Squares (□): Low gas fraction simulations (contain 'fg' in name)

```python
from plotting import plot_depletion_comparison, plot_dtm_dtg_vs_metallicity

# Generate depletion comparison plots using sample data
plot_depletion_comparison('sample_data.h5', output_number=81)

# Generate dust ratio plots using sample data
plot_dtm_dtg_vs_metallicity('sample_data.h5', output_number=81)
```

### Custom Analysis Example

Here's how to read the sample HDF5 file and create custom plots:

```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

def plot_custom_dtm_analysis(h5_filename, output_number=81):
    """
    Custom analysis showing DTM vs initial metallicity with simulation type and gas fraction encoding
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Color mapping for simulation types
    colors = {'G8': 'green', 'G9': 'blue', 'G10': 'purple'}
    
    with h5py.File(h5_filename, 'r') as f:
        for sim_name in f.keys():
            sim_group = f[sim_name]
            
            # Get simulation metadata
            sim_type = sim_group.attrs.get('simulation_type', 'Unknown')
            initial_metallicity = sim_group.attrs.get('initial_metallicity', 1.0)
            
            # Determine marker based on gas fraction
            marker = 's' if 'fg' in sim_name.lower() else 'o'
            gas_fraction_type = 'Low gas fraction' if 'fg' in sim_name.lower() else 'Regular gas fraction'
            
            for depl_config in sim_group.keys():
                depl_group = sim_group[depl_config]
                
                # Get the rmax and zmax values for dust analysis
                rmax_dust = depl_group.attrs.get('rmax_dust_kpc', depl_group.attrs.get('rmax_kpc', 0))
                zmax_dust = depl_group.attrs.get('zmax_dust_kpc', depl_group.attrs.get('zmax_kpc', 0))
                
                # Read dust content data
                dust_key = f'dust_content_{output_number:05d}'
                if dust_key in depl_group:
                    dust_data = depl_group[dust_key][()]
                    
                    # Get column names to find DTM index
                    columns = list(depl_group.attrs['dust_content_columns'])
                    dtm_idx = columns.index('DTM')
                    dtm_value = dust_data[dtm_idx]
                    
                    # Plot with appropriate color and marker
                    ax.scatter(initial_metallicity, dtm_value, 
                             marker=marker, 
                             color=colors.get(sim_type, 'gray'),
                             s=120, alpha=0.7,
                             label=f'{sim_type} - {gas_fraction_type}' if sim_name == list(f.keys())[0] else "")
    
    # Create custom legends
    # Simulation type legend
    sim_legend_elements = [
        plt.scatter([], [], color='green', marker='o', s=120, label='G8'),
        plt.scatter([], [], color='blue', marker='o', s=120, label='G9'),
        plt.scatter([], [], color='purple', marker='o', s=120, label='G10')
    ]
    leg1 = ax.legend(handles=sim_legend_elements, title="Simulation Type", 
                     loc='upper left', frameon=True)
    
    # Gas fraction legend  
    gas_legend_elements = [
        plt.scatter([], [], color='gray', marker='o', s=120, label='Regular gas fraction'),
        plt.scatter([], [], color='gray', marker='s', s=120, label='Low gas fraction (fg)')
    ]
    leg2 = ax.legend(handles=gas_legend_elements, title="Gas Fraction", 
                     loc='upper right', frameon=True)
    ax.add_artist(leg1)  # Add first legend back
    
    ax.set_xlabel('Initial Metallicity (Zsun)')
    ax.set_ylabel('Dust-to-Metal Ratio (DTM)')
    ax.set_title(f'DTM vs Initial Metallicity - Output {output_number:05d}')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'custom_dtm_analysis_output_{output_number:05d}.png', dpi=300, bbox_inches='tight')
    plt.show()

# Usage example with sample data
plot_custom_dtm_analysis('sample_data.h5', output_number=81)
```

### Reading Individual Data Components

```python
# Read specific simulation data from sample dataset
with h5py.File('sample_data.h5', 'r') as f:
    # List all available simulations in the sample dataset
    print("Available simulations:", list(f.keys()))
    
    # Access a specific simulation (note the underscore separator)
    sim_name = list(f.keys())[0]  # Get first available simulation
    sim_group = f[sim_name]
    
    # Get simulation metadata
    sim_type = sim_group.attrs['simulation_type']  # e.g., 'G8'
    initial_metallicity = sim_group.attrs['initial_metallicity']  # e.g., 0.1
    original_path = sim_group.attrs['original_path']
    
    # Access depletion configuration
    depl_config = list(sim_group.keys())[0]  # Get first available configuration
    depl_group = sim_group[depl_config]
    
    # Get region parameters (now separate for depletion and dust)
    rmax_depl = depl_group.attrs['rmax_depl_kpc']
    zmax_depl = depl_group.attrs['zmax_depl_kpc'] 
    rmax_dust = depl_group.attrs['rmax_dust_kpc']
    zmax_dust = depl_group.attrs['zmax_dust_kpc']
    
    # Read depletion vs density data (2D array)
    sample_data = depl_group['output_00081'][()]
    rho = sample_data[:, 0]  # Gas density
    c_depletion = sample_data[:, 1]  # Carbon depletion fraction
    
    # Read integrated dust content data (1D array)
    dust_data = depl_group['dust_content_00081'][()]
    columns = list(depl_group.attrs['dust_content_columns'])
    
    # Extract specific values
    dtm = dust_data[columns.index('DTM')]
    dtg = dust_data[columns.index('DTG')]
    metallicity = dust_data[columns.index('Z')]
    
    print(f"Simulation: {sim_name}")
    print(f"Simulation Type: {sim_type}")
    print(f"Initial Metallicity: {initial_metallicity:.3f} Zsun")
    print(f"DTM: {dtm:.4f}")
    print(f"DTG: {dtg:.6f}")
    print(f"Measured Metallicity: {metallicity:.4f}")
    print(f"Gas fraction type: {'Low' if 'fg' in sim_name.lower() else 'Regular'}")

# Example: Get all simulation types and metallicities from sample data
def get_simulation_summary(h5_filename):
    """Get summary of all simulations in the HDF5 file"""
    with h5py.File(h5_filename, 'r') as f:
        print("Simulation Summary:")
        print("-" * 50)
        for sim_name in f.keys():
            sim_group = f[sim_name]
            sim_type = sim_group.attrs.get('simulation_type', 'Unknown')
            metallicity = sim_group.attrs.get('initial_metallicity', 'Unknown')
            gas_fraction = 'Low' if 'fg' in sim_name.lower() else 'Regular'
            
            print(f"{sim_name}:")
            print(f"  Type: {sim_type}")
            print(f"  Initial Metallicity: {metallicity} Zsun")
            print(f"  Gas Fraction: {gas_fraction}")
            print(f"  Configurations: {list(sim_group.keys())}")
            print()

# Explore the sample dataset
get_simulation_summary('sample_data.h5')
```

### Sample Dataset Contents

The included `sample_data.h5` contains:
- Multiple G8, G9, and G10 simulations
- Various initial metallicity values (0.1, 0.3, 1.0, 3.0 Zsun)
- Both regular and low gas fraction ('fg') models
- Different spatial region configurations
- All processed for output 81 (corresponding to specific simulation time)

### License

This project is licensed under the MIT License. For more details, see the [LICENSE](LICENSE) file.

---

Made with :heart: by [Curro Rodriguez Montero](https://github.com/currodri)

<br>

[Back to top](#top)
