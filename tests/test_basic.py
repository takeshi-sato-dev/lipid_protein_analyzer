#!/usr/bin/env python
"""Basic tests for JOSS submission."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def test_import():
    """Test that the module can be imported."""
    import time_series_analysis
    assert time_series_analysis is not None

def test_constants():
    """Test that required constants exist."""
    import time_series_analysis as tsa
    assert hasattr(tsa, 'TOPOLOGY_FILE')
    assert hasattr(tsa, 'TRAJECTORY_FILE')

def test_functions():
    """Test that main functions exist."""
    import time_series_analysis as tsa
    assert callable(getattr(tsa, 'load_universe'))
    assert callable(getattr(tsa, 'analyze_time_series'))

if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
