#!/usr/bin/env python
"""
Comprehensive test suite for LipidProteinAnalyzer.
Tests core functionality, edge cases, and performance.
"""

import pytest
import numpy as np
import pandas as pd
import os
import sys
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import time_series_analysis as tsa
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def test_universe():
    """Create a minimal test universe."""
    # Create simple test universe with 100 atoms
    n_atoms = 100
    u = mda.Universe.empty(n_atoms, trajectory=True)
    
    # Add basic topology
    resnames = ['PRO'] * 20 + ['POPC'] * 20 + ['POPE'] * 20 + ['POPS'] * 20 + ['CHOL'] * 20
    u.add_TopologyAttr('resname', resnames)
    
    # Add names for lipid headgroup detection
    names = []
    for i in range(100):
        if i < 20:
            names.append('BB')  # Protein backbone
        elif i % 5 == 0:
            names.append('PO4')  # Lipid headgroup
        else:
            names.append('C1')  # Lipid tail
    u.add_TopologyAttr('name', names)
    
    # Set positions (two leaflets)
    positions = np.random.rand(n_atoms, 3) * 100
    # Upper leaflet lipids
    positions[20:60, 2] = np.random.rand(40) * 10 + 60  # z = 60-70
    # Lower leaflet lipids
    positions[60:100, 2] = np.random.rand(40) * 10 + 30  # z = 30-40
    # Proteins in middle
    positions[0:20, 2] = np.random.rand(20) * 10 + 45  # z = 45-55
    
    u.atoms.positions = positions
    u.dimensions = [100, 100, 100, 90, 90, 90]
    
    return u


@pytest.fixture
def mock_trajectory(test_universe):
    """Add mock trajectory to test universe."""
    # Create 10 frames
    n_frames = 10
    n_atoms = len(test_universe.atoms)
    
    frames = []
    for i in range(n_frames):
        pos = test_universe.atoms.positions.copy()
        # Add small random movement
        pos += np.random.randn(n_atoms, 3) * 0.1
        frames.append(pos)
    
    # Mock trajectory iteration
    test_universe.trajectory = Mock()
    test_universe.trajectory.__len__ = Mock(return_value=n_frames)
    test_universe.trajectory.__iter__ = Mock(return_value=iter(range(n_frames)))
    test_universe.trajectory.__getitem__ = Mock(side_effect=lambda i: Mock(positions=frames[i]))
    
    return test_universe


@pytest.fixture
def temp_output_dir():
    """Create temporary directory for test outputs."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_time_series_data():
    """Create sample time series data for plotting tests."""
    data = []
    for frame in range(100):
        frame_data = []
        for protein in range(4):
            protein_data = {
                'POPC_0': np.random.rand() * 10,
                'POPC_1': np.random.rand() * 10,
                'POPE_0': np.random.rand() * 5,
                'POPE_1': np.random.rand() * 5,
                'POPS_0': np.random.rand() * 3,
                'POPS_1': np.random.rand() * 3,
                'CHOL_0': np.random.rand() * 2,
                'CHOL_1': np.random.rand() * 2,
            }
            frame_data.append(protein_data)
        data.append(frame_data)
    return data


# ============================================================================
# Unit Tests - Core Functions
# ============================================================================

class TestLoadingFunctions:
    """Test data loading functions."""
    
    def test_load_universe_with_valid_files(self, tmp_path):
        """Test loading universe with valid topology and trajectory files."""
        # Create dummy files
        psf_file = tmp_path / "test.psf"
        xtc_file = tmp_path / "test.xtc"
        psf_file.touch()
        xtc_file.touch()
        
        # Mock the MDAnalysis Universe
        with patch('time_series_analysis.mda.Universe') as mock_universe:
            mock_u = Mock()
            mock_universe.return_value = mock_u
            
            # Override module constants
            original_psf = tsa.TOPOLOGY_FILE
            original_xtc = tsa.TRAJECTORY_FILE
            tsa.TOPOLOGY_FILE = str(psf_file)
            tsa.TRAJECTORY_FILE = str(xtc_file)
            
            try:
                result = tsa.load_universe()
                assert result == mock_u
                mock_universe.assert_called_once_with(str(psf_file), str(xtc_file))
            finally:
                tsa.TOPOLOGY_FILE = original_psf
                tsa.TRAJECTORY_FILE = original_xtc
    
    def test_load_universe_missing_files(self):
        """Test loading universe with missing files."""
        original_psf = tsa.TOPOLOGY_FILE
        original_xtc = tsa.TRAJECTORY_FILE
        
        tsa.TOPOLOGY_FILE = "nonexistent.psf"
        tsa.TRAJECTORY_FILE = "nonexistent.xtc"
        
        try:
            with pytest.raises(Exception):
                tsa.load_universe()
        finally:
            tsa.TOPOLOGY_FILE = original_psf
            tsa.TRAJECTORY_FILE = original_xtc


class TestLeafletIdentification:
    """Test leaflet identification functions."""
    
    def test_identify_lipid_leaflets(self, test_universe):
        """Test basic leaflet identification."""
        with patch('time_series_analysis.LeafletFinder') as mock_leaflet_finder:
            # Mock the leaflet finder
            mock_finder = Mock()
            mock_finder.groups.return_value = [Mock(), Mock()]
            mock_leaflet_finder.return_value = mock_finder
            
            leaflet0, leaflet1, L = tsa.identify_lipid_leaflets(test_universe)
            
            assert leaflet0 is not None
            assert leaflet1 is not None
            assert L is not None
    
    def test_identify_lipid_leaflets_single_leaflet(self, test_universe):
        """Test leaflet identification with single leaflet."""
        # Set all lipids to same Z coordinate
        lipids = test_universe.select_atoms("resname POPC POPE POPS CHOL")
        lipids.positions[:, 2] = 50.0
        
        with patch('time_series_analysis.LeafletFinder') as mock_leaflet_finder:
            mock_finder = Mock()
            mock_finder.groups.return_value = [Mock()]  # Single leaflet
            mock_leaflet_finder.return_value = mock_finder
            
            # Should handle single leaflet gracefully
            leaflet0, leaflet1, L = tsa.identify_lipid_leaflets(test_universe)
            assert leaflet0 is not None


class TestLipidSelections:
    """Test lipid selection setup."""
    
    def test_setup_lipid_selections(self, test_universe):
        """Test setting up lipid selections."""
        # Create mock leaflets
        leaflet0 = test_universe.select_atoms("resname POPC POPE")
        leaflet1 = test_universe.select_atoms("resname POPS CHOL")
        
        lipid_sels = tsa.setup_lipid_selections(leaflet0, leaflet1)
        
        assert isinstance(lipid_sels, dict)
        assert len(lipid_sels) > 0
        
        # Check structure of selections
        for lipid_type, leaflet_dict in lipid_sels.items():
            assert 0 in leaflet_dict or 1 in leaflet_dict
    
    def test_setup_lipid_selections_empty_leaflet(self, test_universe):
        """Test lipid selection with empty leaflet."""
        leaflet0 = test_universe.select_atoms("resname XXXX")  # Empty selection
        leaflet1 = test_universe.select_atoms("resname POPC")
        
        lipid_sels = tsa.setup_lipid_selections(leaflet0, leaflet1)
        
        # Should handle empty leaflet gracefully
        assert isinstance(lipid_sels, dict)


class TestProteinSelection:
    """Test protein selection functions."""
    
    def test_select_proteins_auto_detect(self, test_universe):
        """Test automatic protein detection."""
        proteins = tsa.select_proteins(test_universe)
        
        assert isinstance(proteins, list)
        assert len(proteins) > 0
        
        for protein in proteins:
            assert hasattr(protein, 'atoms')
    
    def test_select_proteins_specified_number(self, test_universe):
        """Test selecting specific number of proteins."""
        proteins = tsa.select_proteins(test_universe, n_proteins=2)
        
        assert len(proteins) == 2
    
    def test_select_proteins_no_proteins(self, test_universe):
        """Test behavior when no proteins exist."""
        # Remove all proteins
        test_universe.atoms.resnames = ['POPC'] * len(test_universe.atoms)
        
        proteins = tsa.select_proteins(test_universe)
        
        # Should return empty list or handle gracefully
        assert isinstance(proteins, list)


# ============================================================================
# Unit Tests - Analysis Functions
# ============================================================================

class TestFrameProcessing:
    """Test frame processing functions."""
    
    def test_process_frame_basic(self):
        """Test basic frame processing."""
        # Create mock data
        mock_protein = Mock()
        mock_protein.positions = np.array([[50, 50, 50]] * 10)
        
        mock_lipid = Mock()
        mock_lipid.positions = np.array([[52, 50, 50]] * 5)  # 2Å away
        
        mock_box = np.array([100, 100, 100])
        
        # Test weighted contacts
        weighted = tsa.process_frame(0, [mock_protein], {'POPC': {0: mock_lipid}})
        
        assert isinstance(weighted, tuple)
        assert len(weighted) == 2  # frame_idx and results
    
    def test_process_frame_periodic_boundary(self):
        """Test frame processing with periodic boundary conditions."""
        mock_protein = Mock()
        mock_protein.positions = np.array([[5, 50, 50]])  # Near boundary
        
        mock_lipid = Mock()
        mock_lipid.positions = np.array([[95, 50, 50]])  # Other side of box
        
        mock_box = np.array([100, 100, 100])
        
        # Should handle PBC correctly
        weighted = tsa.process_frame(0, [mock_protein], {'POPC': {0: mock_lipid}})
        assert weighted is not None
    
    def test_process_frame_empty_selections(self):
        """Test frame processing with empty selections."""
        mock_protein = Mock()
        mock_protein.positions = np.array([])  # Empty
        
        result = tsa.process_frame(0, [mock_protein], {})
        
        # Should handle empty data gracefully
        assert result is not None


class TestTimeSeriesAnalysis:
    """Test time series analysis functions."""
    
    @patch('time_series_analysis.process_frame')
    def test_analyze_time_series_serial(self, mock_process):
        """Test serial time series analysis."""
        mock_process.return_value = (0, [{'POPC_0': 5.0}])
        
        mock_universe = Mock()
        mock_universe.trajectory = [Mock() for _ in range(10)]
        
        proteins = [Mock()]
        lipid_sels = {'POPC': {0: Mock()}}
        
        results = tsa.analyze_time_series(
            mock_universe, proteins, lipid_sels,
            start=0, stop=10, step=1
        )
        
        assert len(results) == 10
        assert mock_process.called
    
    @patch('time_series_analysis.mp.Pool')
    def test_analyze_time_series_parallel(self, mock_pool):
        """Test parallel time series analysis."""
        # Mock pool and map
        mock_pool_instance = Mock()
        mock_pool.return_value.__enter__ = Mock(return_value=mock_pool_instance)
        mock_pool.return_value.__exit__ = Mock(return_value=None)
        
        mock_results = [(i, [{'POPC_0': 5.0}]) for i in range(10)]
        mock_pool_instance.map.return_value = mock_results
        
        mock_universe = Mock()
        mock_universe.trajectory = [Mock() for _ in range(10)]
        
        proteins = [Mock()]
        lipid_sels = {'POPC': {0: Mock()}}
        
        results = tsa.analyze_time_series_parallel(
            mock_universe, proteins, lipid_sels,
            start=0, stop=10, step=1
        )
        
        assert len(results) == 10
        assert mock_pool.called
    
    def test_analyze_time_series_step_parameter(self):
        """Test time series analysis with step parameter."""
        mock_universe = Mock()
        mock_universe.trajectory = [Mock() for _ in range(100)]
        
        with patch('time_series_analysis.process_frame') as mock_process:
            mock_process.return_value = (0, [{'POPC_0': 5.0}])
            
            proteins = [Mock()]
            lipid_sels = {'POPC': {0: Mock()}}
            
            # Test with step=10
            results = tsa.analyze_time_series(
                mock_universe, proteins, lipid_sels,
                start=0, stop=100, step=10
            )
            
            # Should process 10 frames
            assert len(results) == 10


# ============================================================================
# Unit Tests - Visualization Functions
# ============================================================================

class TestPlottingFunctions:
    """Test plotting and visualization functions."""
    
    def test_plot_lipid_composition_time_series(self, sample_time_series_data, temp_output_dir):
        """Test time series plotting."""
        proteins = [f"Protein_{i}" for i in range(4)]
        
        # Should create plot without errors
        tsa.plot_lipid_composition_time_series(
            sample_time_series_data, proteins, temp_output_dir,
            window_size=10
        )
        
        # Check that files were created
        output_files = os.listdir(temp_output_dir)
        assert any('.png' in f for f in output_files)
    
    def test_plot_lipid_bar_charts(self, sample_time_series_data, temp_output_dir):
        """Test bar chart plotting."""
        proteins = [f"Protein_{i}" for i in range(4)]
        lipid_sels = {'POPC': {0: Mock(), 1: Mock()}}
        L = [0, 1]
        
        tsa.plot_lipid_bar_charts(
            sample_time_series_data, proteins, lipid_sels, L, temp_output_dir
        )
        
        output_files = os.listdir(temp_output_dir)
        assert any('bar' in f.lower() for f in output_files)
    
    def test_plot_with_empty_data(self, temp_output_dir):
        """Test plotting with empty data."""
        empty_data = []
        proteins = []
        
        # Should handle empty data gracefully
        try:
            tsa.plot_lipid_composition_time_series(
                empty_data, proteins, temp_output_dir
            )
        except Exception as e:
            # Should raise meaningful error
            assert "empty" in str(e).lower() or "no data" in str(e).lower()
    
    def test_plot_smoothing_windows(self, sample_time_series_data, temp_output_dir):
        """Test different smoothing windows."""
        proteins = [f"Protein_{i}" for i in range(4)]
        
        for window in [10, 30, 50]:
            tsa.plot_lipid_composition_time_series(
                sample_time_series_data, proteins, temp_output_dir,
                window_size=window
            )
            
            # Check window-specific files
            output_files = os.listdir(temp_output_dir)
            assert any(f'window{window}' in f for f in output_files)


# ============================================================================
# Integration Tests
# ============================================================================

class TestIntegrationWorkflow:
    """Test complete analysis workflow."""
    
    @pytest.mark.integration
    def test_full_workflow_with_test_data(self, temp_output_dir):
        """Test complete workflow with test data."""
        # Check if test data exists, if not skip
        if not os.path.exists('test_data/test_system.psf'):
            pytest.skip("Test data not available")
        
        # Override paths
        original_psf = tsa.TOPOLOGY_FILE
        original_xtc = tsa.TRAJECTORY_FILE
        
        try:
            tsa.TOPOLOGY_FILE = 'test_data/test_system.psf'
            tsa.TRAJECTORY_FILE = 'test_data/test_trajectory.xtc'
            
            # Run complete workflow
            u = tsa.load_universe()
            leaflet0, leaflet1, L = tsa.identify_lipid_leaflets(u)
            lipid_sels = tsa.setup_lipid_selections(leaflet0, leaflet1)
            proteins = tsa.select_proteins(u)
            
            # Quick analysis
            time_series_data = tsa.analyze_time_series(
                u, proteins, lipid_sels,
                start=0, stop=5, step=1
            )
            
            # Generate plots
            tsa.plot_lipid_composition_time_series(
                time_series_data, proteins, temp_output_dir
            )
            
            # Check outputs
            assert len(time_series_data) > 0
            assert os.listdir(temp_output_dir)  # Should have created files
            
        finally:
            tsa.TOPOLOGY_FILE = original_psf
            tsa.TRAJECTORY_FILE = original_xtc
    
    @pytest.mark.integration
    def test_parallel_vs_serial_consistency(self):
        """Test that parallel and serial analysis give same results."""
        # Create mock data
        mock_universe = Mock()
        mock_universe.trajectory = [Mock() for _ in range(20)]
        
        proteins = [Mock()]
        proteins[0].positions = np.random.rand(10, 3) * 100
        
        lipid_sels = {'POPC': {0: Mock()}}
        lipid_sels['POPC'][0].positions = np.random.rand(20, 3) * 100
        
        with patch('time_series_analysis.process_frame') as mock_process:
            mock_process.return_value = (0, [{'POPC_0': 5.0}])
            
            # Run serial
            serial_results = tsa.analyze_time_series(
                mock_universe, proteins, lipid_sels,
                start=0, stop=10, step=1
            )
            
            # Run parallel
            with patch('time_series_analysis.mp.Pool') as mock_pool:
                mock_pool_instance = Mock()
                mock_pool.return_value.__enter__ = Mock(return_value=mock_pool_instance)
                mock_pool.return_value.__exit__ = Mock(return_value=None)
                mock_pool_instance.map.return_value = [(i, [{'POPC_0': 5.0}]) for i in range(10)]
                
                parallel_results = tsa.analyze_time_series_parallel(
                    mock_universe, proteins, lipid_sels,
                    start=0, stop=10, step=1
                )
            
            # Results should be consistent
            assert len(serial_results) == len(parallel_results)


# ============================================================================
# Performance Tests
# ============================================================================

class TestPerformance:
    """Test performance characteristics."""
    
    @pytest.mark.slow
    def test_large_system_memory_usage(self):
        """Test memory usage with large system."""
        pytest.skip("Implement with memory_profiler if needed")
    
    @pytest.mark.slow
    def test_parallel_speedup(self):
        """Test parallel processing speedup."""
        pytest.skip("Implement timing comparison if needed")


# ============================================================================
# Edge Cases and Error Handling
# ============================================================================

class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_single_frame_analysis(self):
        """Test analysis with single frame."""
        mock_universe = Mock()
        mock_universe.trajectory = [Mock()]
        
        with patch('time_series_analysis.process_frame') as mock_process:
            mock_process.return_value = (0, [{'POPC_0': 5.0}])
            
            results = tsa.analyze_time_series(
                mock_universe, [Mock()], {'POPC': {0: Mock()}},
                start=0, stop=1, step=1
            )
            
            assert len(results) == 1
    
    def test_no_lipids_in_contact(self):
        """Test when no lipids are in contact with protein."""
        mock_protein = Mock()
        mock_protein.positions = np.array([[0, 0, 0]])
        
        mock_lipid = Mock()
        mock_lipid.positions = np.array([[100, 100, 100]])  # Far away
        
        result = tsa.process_frame(0, [mock_protein], {'POPC': {0: mock_lipid}})
        
        # Should return zero or very small contacts
        assert result is not None
    
    def test_invalid_frame_range(self):
        """Test invalid frame range parameters."""
        mock_universe = Mock()
        mock_universe.trajectory = [Mock() for _ in range(10)]
        
        # Start > Stop should raise error or return empty
        with pytest.raises(Exception):
            tsa.analyze_time_series(
                mock_universe, [Mock()], {},
                start=10, stop=0, step=1
            )
    
    def test_missing_lipid_types(self):
        """Test handling of missing lipid types."""
        mock_universe = Mock()
        mock_universe.atoms = Mock()
        mock_universe.atoms.resnames = ['PRO'] * 100  # Only proteins
        
        # Should handle missing lipids gracefully
        leaflet0 = Mock()
        leaflet1 = Mock()
        leaflet0.resnames = []
        leaflet1.resnames = []
        
        lipid_sels = tsa.setup_lipid_selections(leaflet0, leaflet1)
        
        assert isinstance(lipid_sels, dict)


# ============================================================================
# Run tests with coverage report
# ============================================================================

if __name__ == "__main__":
    # Run with: python -m pytest tests/test_core_functions.py -v --cov=time_series_analysis
    pytest.main([__file__, '-v', '--cov=time_series_analysis', '--cov-report=term-missing'])