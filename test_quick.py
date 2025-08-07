#!/usr/bin/env python
"""
Quick test to verify the analysis works with test data.
Run this after installation to check everything works.

Usage:
    python test_quick.py
"""

import sys
import os
import time

def run_quick_test():
    """Run a quick test with the test data."""
    
    print("="*60)
    print("LipidProteinAnalyzer - Quick Test")
    print("="*60)
    
    # Import the analysis module
    try:
        import time_series_analysis as tsa
        print("✓ Module imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import module: {e}")
        print("\nPlease install the package:")
        print("  pip install -e .")
        return False
    
    # Check test data exists
    if not os.path.exists('test_data/test_system.psf'):
        print("\n✗ Test data not found!")
        print("Please run: python create_test_data.py")
        return False
    
    print("✓ Test data found")
    
    # Override configuration for test data
    print("\nConfiguring for test data...")
    tsa.TOPOLOGY_FILE = 'test_data/test_system.psf'
    tsa.TRAJECTORY_FILE = 'test_data/test_trajectory.xtc'
    tsa.OUTPUT_DIR = 'test_output'
    
    # Clean up old output
    os.makedirs('test_output', exist_ok=True)
    
    start_time = time.time()
    
    try:
        # Step 1: Load data
        print("\n1. Loading universe...")
        u = tsa.load_universe()
        print(f"   Loaded: {len(u.atoms)} atoms, {len(u.trajectory)} frames")
        
        # Step 2: Identify leaflets
        print("\n2. Identifying leaflets...")
        leaflet0, leaflet1, L = tsa.identify_lipid_leaflets(u)
        
        # Step 3: Setup selections
        print("\n3. Setting up selections...")
        lipid_sels = tsa.setup_lipid_selections(leaflet0, leaflet1)
        proteins = tsa.select_proteins(u)  # Auto-detect proteins (no n_proteins argument)
        print(f"   Found: {len(proteins)} proteins, {len(lipid_sels)} lipid types")
        
        # Step 4: Run quick analysis (10 frames only for speed)
        print("\n4. Running analysis (10 frames)...")
        time_series_data = tsa.analyze_time_series(
            u, proteins, lipid_sels,
            start=0, stop=10, step=1, interval=1
        )
        print(f"   Analyzed: {len(time_series_data)} frames")
        
        # Step 5: Generate one plot
        print("\n5. Generating sample plot...")
        tsa.plot_lipid_composition_time_series(
            time_series_data, proteins, 'test_output', window_size=3
        )
        
        # Check output
        output_files = [f for f in os.listdir('test_output') if f.endswith('.png')]
        print(f"   Generated: {len(output_files)} plots")
        
        elapsed = time.time() - start_time
        
        print("\n" + "="*60)
        print(f"✓ TEST PASSED in {elapsed:.1f} seconds!")
        print("="*60)
        print("\nOutput files in: test_output/")
        print("Sample files generated:")
        for f in output_files[:3]:  # Show first 3 files
            print(f"  - {f}")
        
        return True
        
    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        print("\nDebugging info:")
        print(f"  Current directory: {os.getcwd()}")
        print(f"  Python version: {sys.version}")
        
        # Try to provide helpful error messages
        if "No module named" in str(e):
            print("\n  Missing dependency. Install with:")
            print("    pip install -r requirements.txt")
        elif "File" in str(e) and "not found" in str(e):
            print("\n  Data files missing. Run:")
            print("    python create_test_data.py")
        
        return False

if __name__ == "__main__":
    success = run_quick_test()
    sys.exit(0 if success else 1)