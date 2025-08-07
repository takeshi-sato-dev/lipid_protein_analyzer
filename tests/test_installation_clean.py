#!/usr/bin/env python
"""
Clean installation tests without pytest warnings.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def test_import():
    """Test that the module can be imported."""
    import time_series_analysis
    assert time_series_analysis is not None
    # No return statement needed

def test_functions():
    """Test that main functions exist."""
    import time_series_analysis as tsa
    
    required_functions = [
        'load_universe',
        'identify_lipid_leaflets', 
        'setup_lipid_selections',
        'select_proteins',
        'analyze_time_series',
        'analyze_time_series_parallel'
    ]
    
    for func in required_functions:
        assert hasattr(tsa, func), f"Missing function: {func}"
        assert callable(getattr(tsa, func)), f"Not callable: {func}"
    # No return statement needed

def test_constants():
    """Test that required constants are defined."""
    import time_series_analysis as tsa
    
    required_constants = [
        'TOPOLOGY_FILE',
        'TRAJECTORY_FILE',
        'OUTPUT_DIR',
        'START',
        'STOP', 
        'STEP',
        'CONTACT_CUTOFF'
    ]
    
    for const in required_constants:
        assert hasattr(tsa, const), f"Missing constant: {const}"
    # No return statement needed

def test_dependencies():
    """Test that required dependencies can be imported."""
    try:
        import MDAnalysis
        import numpy
        import pandas
        import matplotlib
        
        assert MDAnalysis is not None
        assert numpy is not None
        assert pandas is not None
        assert matplotlib is not None
    except ImportError as e:
        assert False, f"Missing dependency: {e}"
    # No return statement needed

if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
