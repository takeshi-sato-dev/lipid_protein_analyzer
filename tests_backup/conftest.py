"""
Pytest configuration and shared fixtures for LipidProteinAnalyzer tests.
"""

import os
import sys
import pytest
import tempfile
import shutil
from pathlib import Path

# Add parent directory to Python path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ============================================================================
# Configuration
# ============================================================================

def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "visualization: marks tests that generate plots"
    )
    config.addinivalue_line(
        "markers", "requires_data: marks tests that require test data files"
    )


# ============================================================================
# Session Fixtures
# ============================================================================

@pytest.fixture(scope="session")
def test_data_dir():
    """Path to test data directory."""
    return Path(__file__).parent / "test_data"


@pytest.fixture(scope="session", autouse=True)
def setup_test_data(test_data_dir):
    """Ensure test data exists before running tests."""
    if not test_data_dir.exists():
        test_data_dir.mkdir(parents=True, exist_ok=True)
        
        # Try to generate test data
        try:
            from test_data_generator import main as generate_data
            print("\nGenerating test data...")
            result = generate_data()
            if result != 0:
                pytest.skip("Could not generate test data")
        except ImportError:
            print("Warning: test_data_generator.py not found")
            # Continue anyway - some tests may not need data
    
    yield
    
    # Cleanup is optional - preserve test data between runs


@pytest.fixture(scope="session")
def test_psf_file(test_data_dir):
    """Path to test PSF file."""
    psf_file = test_data_dir / "test_system.psf"
    if not psf_file.exists():
        pytest.skip(f"Test PSF file not found: {psf_file}")
    return str(psf_file)


@pytest.fixture(scope="session")
def test_xtc_file(test_data_dir):
    """Path to test XTC file."""
    xtc_file = test_data_dir / "test_trajectory.xtc"
    if not xtc_file.exists():
        pytest.skip(f"Test XTC file not found: {xtc_file}")
    return str(xtc_file)


# ============================================================================
# Function Fixtures
# ============================================================================

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    temp_path = tempfile.mkdtemp(prefix="lpa_test_")
    yield Path(temp_path)
    # Cleanup
    shutil.rmtree(temp_path, ignore_errors=True)


@pytest.fixture
def mock_time_series_data():
    """Generate mock time series data for testing."""
    import numpy as np
    
    n_frames = 100
    n_proteins = 4
    lipid_types = ['POPC', 'POPE', 'POPS', 'CHOL']
    
    data = []
    for frame in range(n_frames):
        frame_data = []
        for protein in range(n_proteins):
            protein_data = {}
            for lipid in lipid_types:
                # Random counts for each leaflet
                protein_data[f'{lipid}_0'] = np.random.poisson(5)
                protein_data[f'{lipid}_1'] = np.random.poisson(5)
            frame_data.append(protein_data)
        data.append(frame_data)
    
    return data


@pytest.fixture
def analysis_module():
    """Import and return the analysis module."""
    import time_series_analysis
    return time_series_analysis


@pytest.fixture
def mock_universe():
    """Create a mock MDAnalysis Universe."""
    import MDAnalysis as mda
    import numpy as np
    
    # Create simple universe
    n_atoms = 100
    u = mda.Universe.empty(n_atoms, trajectory=True)
    
    # Add topology
    resnames = ['PRO'] * 20 + ['POPC'] * 20 + ['POPE'] * 20 + ['POPS'] * 20 + ['CHOL'] * 20
    u.add_TopologyAttr('resname', resnames)
    
    names = []
    for i in range(100):
        if i < 20:
            names.append('BB')
        elif i % 5 == 0:
            names.append('PO4')
        else:
            names.append('C1')
    u.add_TopologyAttr('name', names)
    
    # Set positions
    positions = np.random.rand(n_atoms, 3) * 100
    u.atoms.positions = positions
    u.dimensions = [100, 100, 100, 90, 90, 90]
    
    return u


# ============================================================================
# Fixture for Module Path Override
# ============================================================================

@pytest.fixture
def override_paths(analysis_module, test_psf_file, test_xtc_file):
    """Override module paths for testing."""
    # Save original values
    original_psf = analysis_module.TOPOLOGY_FILE
    original_xtc = analysis_module.TRAJECTORY_FILE
    original_output = analysis_module.OUTPUT_DIR
    
    # Override with test values
    analysis_module.TOPOLOGY_FILE = test_psf_file
    analysis_module.TRAJECTORY_FILE = test_xtc_file
    analysis_module.OUTPUT_DIR = "test_output"
    
    yield analysis_module
    
    # Restore original values
    analysis_module.TOPOLOGY_FILE = original_psf
    analysis_module.TRAJECTORY_FILE = original_xtc
    analysis_module.OUTPUT_DIR = original_output


# ============================================================================
# Performance Testing Fixtures
# ============================================================================

@pytest.fixture
def benchmark_data():
    """Generate larger dataset for benchmarking."""
    import numpy as np
    
    # Larger dataset for performance testing
    n_frames = 1000
    n_proteins = 4
    n_lipids = 200
    
    data = {
        'n_frames': n_frames,
        'n_proteins': n_proteins,
        'n_lipids': n_lipids,
        'positions': np.random.rand(n_frames, n_lipids, 3) * 100
    }
    
    return data


# ============================================================================
# Utility Functions
# ============================================================================

def pytest_collection_modifyitems(config, items):
    """Modify test collection to add markers based on test names."""
    for item in items:
        # Auto-mark slow tests
        if "slow" in item.nodeid:
            item.add_marker(pytest.mark.slow)
        
        # Auto-mark visualization tests
        if "plot" in item.nodeid or "visualiz" in item.nodeid:
            item.add_marker(pytest.mark.visualization)
        
        # Auto-mark integration tests
        if "integration" in item.nodeid or "workflow" in item.nodeid:
            item.add_marker(pytest.mark.integration)


@pytest.fixture(autouse=True)
def reset_matplotlib():
    """Reset matplotlib to prevent figure accumulation."""
    import matplotlib.pyplot as plt
    yield
    plt.close('all')


# ============================================================================
# Command Line Options
# ============================================================================

def pytest_addoption(parser):
    """Add custom command line options."""
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="run slow tests"
    )
    parser.addoption(
        "--generate-plots",
        action="store_true",
        default=False,
        help="generate and save plots during tests"
    )


def pytest_runtest_setup(item):
    """Skip tests based on command line options."""
    if 'slow' in item.keywords and not item.config.getoption("--run-slow"):
        pytest.skip("need --run-slow option to run slow tests")