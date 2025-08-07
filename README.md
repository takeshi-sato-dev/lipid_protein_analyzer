# LipidProteinAnalyzer

[![DOI](https://doi.org/10.5281/zenodo.16756837)

Time-series analysis tool for lipid-protein interactions in molecular dynamics simulations with parallel processing support.

## Features

- Time-resolved analysis of lipid-protein contacts
- Parallel processing for large trajectories
- Distance-weighted and simple counting metrics
- Automated leaflet identification
- Multiple visualization options with smoothing
- Support for Martini coarse-grained simulations

## Installation

```bash
git clone https://github.com/yourusername/lipid_protein_analyzer.git
cd lipid_protein_analyzer
pip install -r requirements.txt
pip install .
```

## Usage

### As a Python script

```python
# Run directly
python time_series_analysis.py

# Or import functions
from time_series_analysis import load_universe, analyze_time_series_parallel
from time_series_analysis import plot_lipid_composition_time_series

# Load trajectory
u = load_universe()

# Setup system
leaflet0, leaflet1, L = identify_lipid_leaflets(u)
lipid_sels = setup_lipid_selections(leaflet0, leaflet1)
proteins = select_proteins(u, n_proteins=4)

# Run analysis
time_series_data = analyze_time_series_parallel(
    u, proteins, lipid_sels, 
    start=60000, stop=80000, step=50
)

# Generate plots
plot_lipid_composition_time_series(
    time_series_data, proteins, 
    output_dir='results/', window_size=20
)
```

### Configuration

Modify the constants at the top of `time_series_analysis.py`:

```python
START = 60000  # Starting frame
STOP = 80000   # Ending frame
STEP = 50      # Frame step
CONTACT_CUTOFF = 8.0  # Angstroms
TOPOLOGY_FILE = 'step5_assembly.psf'
TRAJECTORY_FILE = 'md_wrapped.xtc'
OUTPUT_DIR = 'lipid_protein_timeseries'
```

## System Requirements

- Python 3.8+
- 4+ GB RAM for large trajectories
- Multiple CPU cores recommended for parallel processing

## Output

The analysis generates:
- Time series plots with smoothing options (10, 30, 50 frame windows)
- Bar charts of average lipid composition
- Both distance-weighted and simple count metrics
- Leaflet-specific analysis
- Output in PNG and PS formats

## Citation

If you use this software, please cite:

```bibtex
@article{YourName2025,
  title={LipidProteinAnalyzer: Time-series analysis of lipid-protein interactions in molecular dynamics simulations},
  author={Takeshi Sato},
  journal={Journal of Open Source Software},
  year={2025},
  doi={10.21105/joss.XXXXX}
}
```

## License

MIT License - see [LICENSE](LICENSE) file