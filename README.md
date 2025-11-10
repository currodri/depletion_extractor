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
├── simulation_1/
│   └── DEPLETION_rmax{X}kpc_zmax{Y}kpc/
│       ├── output_{NNNNN}           # 2D array: depletion data vs density
│       ├── dust_content_{NNNNN}     # 1D array: integrated dust/gas properties
│       └── attributes:
│           ├── rmax_kpc             # Radial extent in kpc
│           ├── zmax_kpc             # Vertical extent in kpc  
│           └── dust_content_columns # Column names for dust_content data
├── simulation_2/
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

### Input File Format

Create a text file listing simulation paths with their corresponding rmax and zmax values:

```
# path_list.txt
/path/to/simulation1 4.0 1.0
/path/to/simulation2 8.0 0.2  
/path/to/simulation3 12.0 0.5
```

Each line contains: `simulation_path rmax_value zmax_value`

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
# -f output_data.h5 : Output HDF5 filename (optional, default: depletion_data.h5)
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

## :chart_with_upwards_trend: Example: Reading and Plotting DTM Data ##

Here's how to read the extracted HDF5 file and plot the dust-to-metal ratios for all simulations:

```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

def plot_dtm_comparison(h5_filename, output_number=81):
    """
    Plot dust-to-metal ratios for all simulations and depletion configurations
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    with h5py.File(h5_filename, 'r') as f:
        for sim_name in f.keys():
            sim_group = f[sim_name]
            
            for depl_config in sim_group.keys():
                depl_group = sim_group[depl_config]
                
                # Get the rmax and zmax values
                rmax = depl_group.attrs['rmax_kpc']
                zmax = depl_group.attrs['zmax_kpc']
                
                # Read dust content data
                dust_key = f'dust_content_{output_number:05d}'
                if dust_key in depl_group:
                    dust_data = depl_group[dust_key][()]
                    
                    # Get column names to find DTM index
                    columns = list(depl_group.attrs['dust_content_columns'])
                    dtm_idx = columns.index('DTM')
                    dtm_value = dust_data[dtm_idx]
                    
                    # Create label and plot
                    label = f'{sim_name} (r={rmax:.1f}, z={zmax:.1f})'
                    ax.scatter(rmax, dtm_value, label=label, s=100, alpha=0.7)
    
    ax.set_xlabel('Radial Extent (rmax) [kpc]')
    ax.set_ylabel('Dust-to-Metal Ratio (DTM)')
    ax.set_title(f'Dust-to-Metal Ratios - Output {output_number:05d}')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'dtm_comparison_output_{output_number:05d}.png', dpi=300, bbox_inches='tight')
    plt.show()

# Usage example
plot_dtm_comparison('output_data.h5', output_number=81)
```

### Reading Individual Data Components

```python
# Read specific simulation data
with h5py.File('output_data.h5', 'r') as f:
    # Access a specific simulation and depletion configuration
    sim_group = f['simulation_name']
    depl_group = sim_group['DEPLETION_rmax4.0kpc_zmax1.0kpc']
    
    # Read depletion vs density data (2D array)
    depletion_data = depl_group['output_00081'][()]
    rho = depletion_data[:, 0]  # Gas density
    c_depletion = depletion_data[:, 1]  # Carbon depletion fraction
    
    # Read integrated dust content data (1D array)
    dust_data = depl_group['dust_content_00081'][()]
    columns = list(depl_group.attrs['dust_content_columns'])
    
    # Extract specific values
    dtm = dust_data[columns.index('DTM')]
    dtg = dust_data[columns.index('DTG')]
    metallicity = dust_data[columns.index('Z')]
    
    print(f"DTM: {dtm:.4f}")
    print(f"DTG: {dtg:.6f}")
    print(f"Metallicity: {metallicity:.4f}")
```

## :memo: License ##

This project is under license from MIT. For more details, see the [LICENSE](LICENSE) file.


Made with :heart: by <a href="https://github.com/currodri" target="_blank">Curro Rodriguez Montero</a>

&#xa0;

<a href="#top">Back to top</a>
