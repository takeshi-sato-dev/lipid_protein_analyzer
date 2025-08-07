---
title: 'LipidProteinAnalyzer: Parallel time-series analysis of lipid-protein interactions in molecular dynamics simulations'
tags:
  - Python
  - molecular dynamics
  - lipid-protein interactions
  - membrane proteins
  - parallel processing
authors:
  - name: Takeshi Sato
    orcid: 0009-0006-9156-8655
    affiliation: 1
affiliations:
 - name: Kyoto Pharmaceutical University
   index: 1
date: 7 August 2025
bibliography: paper.bib
---

# Summary

`LipidProteinAnalyzer` is a Python tool for analyzing time-dependent lipid-protein interactions in molecular dynamics (MD) simulations. The software automatically identifies membrane leaflets, calculates distance-weighted and binary contact metrics, and generates publication-ready visualizations with customizable smoothing windows. It features parallel processing capabilities that reduce analysis time from hours to minutes for large trajectories, making it practical for analyzing microsecond-scale simulations.

# Statement of need

Membrane proteins function within specific lipid environments where individual lipid species can modulate protein structure and function. While general-purpose MD analysis tools exist [@Michaud-Agrawal2011], analyzing time-resolved lipid-protein interactions typically requires custom scripts that are difficult to reproduce and optimize. Researchers need tools that can efficiently process large trajectories (>100 GB) while tracking how lipid compositions change around proteins over time.

Current challenges include:
- Processing large trajectories is computationally expensive without parallelization
- Comparing distance-weighted versus simple counting metrics requires separate implementations
- Generating consistent visualizations across different systems and time windows
- Distinguishing upper and lower leaflet interactions in membrane systems

`LipidProteinAnalyzer` addresses these needs with:
- **Parallel processing** using Python's multiprocessing, achieving near-linear speedup with up to 75% CPU utilization
- **Dual metrics**: both distance-weighted contacts (w = e^(-d/2)) and simple binary counting
- **Automated visualization** with 10, 30, and 50-frame smoothing windows
- **Leaflet-specific analysis** using MDAnalysis's LeafletFinder algorithm
- **Memory-efficient processing** through frame chunking


# Implementation and Performance

The tool processes trajectories in parallel by distributing frames across multiple CPU cores. For a typical system (4 proteins, ~50,000 atoms, 20,000 frames), performance benchmarks show:
- Serial processing: ~2.5 frames/second
- Parallel (8 cores): ~18 frames/second
- Memory usage: <2 GB regardless of trajectory size

Key features include:
- Vectorized NumPy operations for distance calculations with periodic boundary conditions
- Automatic detection and use of 75% of available CPU cores
- Both weighted (distance-dependent) and simple (binary) contact metrics
- Multiple visualization outputs with different smoothing windows for publication flexibility

![Figure 1: Time-series analysis of lipid-protein interactions. The software tracks multiple lipid types (POPC, POPE, POPS, and cholesterol) around membrane proteins over the simulation trajectory. Each panel shows one protein with smoothed trends highlighting temporal changes in lipid composition.](figures/figure1.png)

The figure above demonstrates the visualization capabilities of LipidProteinAnalyzer, showing how different lipid species interact with membrane proteins over time. The automatic smoothing feature helps identify trends while preserving the underlying data structure.

# Usage

Basic usage with default parameters:

```python
from time_series_analysis import *

# Load system
u = load_universe()
leaflet0, leaflet1, L = identify_lipid_leaflets(u)
lipid_sels = setup_lipid_selections(leaflet0, leaflet1)
proteins = select_proteins(u, n_proteins=4)

# Run parallel analysis
time_series_data = analyze_time_series_parallel(
    u, proteins, lipid_sels, 
    start=60000, stop=80000, step=50
)

# Generate visualizations
plot_lipid_composition_time_series(time_series_data, proteins, 'output/')
plot_lipid_bar_charts(time_series_data, proteins, lipid_sels, L, 'output/')
```

The software has been tested on Martini coarse-grained simulations but supports any trajectory format compatible with MDAnalysis.

# Acknowledgements

We acknowledge the MDAnalysis community for their trajectory analysis framework.

# References