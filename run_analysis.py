#!/usr/bin/env python
"""
Wrapper script to run time_series_analysis with custom parameters.
This makes it easier for users to run the analysis without modifying the main code.
"""

import argparse
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description='Run lipid-protein interaction time series analysis'
    )
    
    parser.add_argument(
        '--start', type=int, default=60000,
        help='Starting frame (default: 60000)'
    )
    parser.add_argument(
        '--stop', type=int, default=80000,
        help='Ending frame (default: 80000)'
    )
    parser.add_argument(
        '--step', type=int, default=50,
        help='Frame step (default: 50)'
    )
    parser.add_argument(
        '--topology', type=str, default='step5_assembly.psf',
        help='Topology file (default: step5_assembly.psf)'
    )
    parser.add_argument(
        '--trajectory', type=str, default='md_wrapped.xtc',
        help='Trajectory file (default: md_wrapped.xtc)'
    )
    parser.add_argument(
        '--output', type=str, default='lipid_protein_timeseries',
        help='Output directory (default: lipid_protein_timeseries)'
    )
    parser.add_argument(
        '--parallel', action='store_true', default=True,
        help='Use parallel processing (default: True)'
    )
    parser.add_argument(
        '--no-parallel', dest='parallel', action='store_false',
        help='Disable parallel processing'
    )
    
    args = parser.parse_args()
    
    # Import the analysis functions
    try:
        from time_series_analysis import (
            load_universe,
            identify_lipid_leaflets,
            setup_lipid_selections,
            select_proteins,
            analyze_time_series_parallel,
            analyze_time_series,
            plot_lipid_composition_time_series,
            plot_lipid_bar_charts,
            plot_simple_time_series,
            plot_simple_bar_charts
        )
    except ImportError as e:
        print(f"Error: Could not import time_series_analysis: {e}")
        print("Make sure the package is installed: pip install -e .")
        return 1
    
    # Override the constants (monkey patching for this run only)
    import time_series_analysis as tsa
    tsa.START = args.start
    tsa.STOP = args.stop
    tsa.STEP = args.step
    tsa.TOPOLOGY_FILE = args.topology
    tsa.TRAJECTORY_FILE = args.trajectory
    tsa.OUTPUT_DIR = args.output
    
    print(f"Running analysis with parameters:")
    print(f"  Frames: {args.start} to {args.stop} (step {args.step})")
    print(f"  Topology: {args.topology}")
    print(f"  Trajectory: {args.trajectory}")
    print(f"  Output: {args.output}")
    print(f"  Parallel: {args.parallel}")
    
    try:
        # Run the analysis
        u = load_universe()
        leaflet0, leaflet1, L = identify_lipid_leaflets(u)
        lipid_sels = setup_lipid_selections(leaflet0, leaflet1)
        proteins = select_proteins(u)
        
        if args.parallel:
            time_series_data = analyze_time_series_parallel(
                u, proteins, lipid_sels, 
                args.start, args.stop, args.step
            )
        else:
            time_series_data = analyze_time_series(
                u, proteins, lipid_sels,
                args.start, args.stop, args.step
            )
        
        # Generate plots
        output_dir = Path(args.output)
        output_dir.mkdir(exist_ok=True)
        
        print("Generating visualizations...")
        for window_size in [10, 30, 50]:
            plot_lipid_composition_time_series(
                time_series_data, proteins, output_dir, 
                window_size=window_size
            )
            plot_simple_time_series(
                time_series_data, proteins, output_dir,
                window_size=window_size
            )
        
        plot_lipid_bar_charts(time_series_data, proteins, lipid_sels, L, output_dir)
        plot_simple_bar_charts(time_series_data, proteins, lipid_sels, L, output_dir)
        
        print(f"Analysis complete! Results saved to {output_dir}")
        return 0
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())