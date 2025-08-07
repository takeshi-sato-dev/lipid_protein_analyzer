#!/usr/bin/env python
"""
Simple test script to verify installation for JOSS submission.
Run this after installing the package to ensure everything works.
"""

import sys

def test_import():
    """Test that the module can be imported."""
    try:
        import time_series_analysis
        print("✓ Module import successful")
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False

def test_functions():
    """Test that main functions are available."""
    try:
        from time_series_analysis import (
            load_universe,
            identify_lipid_leaflets,
            setup_lipid_selections,
            select_proteins,
            analyze_time_series_parallel,
            plot_lipid_composition_time_series
        )
        print("✓ All main functions are importable")
        return True
    except ImportError as e:
        print(f"✗ Function import failed: {e}")
        return False

def test_constants():
    """Test that configuration constants are defined."""
    try:
        from time_series_analysis import (
            START, STOP, STEP,
            CONTACT_CUTOFF,
            TOPOLOGY_FILE,
            TRAJECTORY_FILE,
            OUTPUT_DIR
        )
        print("✓ Configuration constants are defined")
        print(f"  - Analysis range: {START} to {STOP} with step {STEP}")
        print(f"  - Contact cutoff: {CONTACT_CUTOFF} Å")
        print(f"  - Output directory: {OUTPUT_DIR}")
        return True
    except ImportError as e:
        print(f"✗ Constants not found: {e}")
        return False

def test_dependencies():
    """Test that required dependencies are installed."""
    dependencies = [
        'MDAnalysis',
        'numpy',
        'pandas',
        'matplotlib',
        'seaborn',
        'tqdm'
    ]
    
    all_ok = True
    for dep in dependencies:
        try:
            __import__(dep)
            print(f"✓ {dep} is installed")
        except ImportError:
            print(f"✗ {dep} is NOT installed")
            all_ok = False
    
    return all_ok

def main():
    """Run all tests."""
    print("=" * 50)
    print("LipidProteinAnalyzer Installation Test")
    print("=" * 50)
    
    tests = [
        ("Import test", test_import),
        ("Function availability", test_functions),
        ("Configuration constants", test_constants),
        ("Dependencies", test_dependencies)
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        results.append(test_func())
    
    print("\n" + "=" * 50)
    if all(results):
        print("✓ All tests passed! Package is ready for JOSS submission.")
        return 0
    else:
        print("✗ Some tests failed. Please fix the issues above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())