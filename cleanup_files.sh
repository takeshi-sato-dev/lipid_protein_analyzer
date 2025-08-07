#!/bin/bash
# cleanup_files.sh - Fix file names and structure

echo "Cleaning up file structure..."

# Create tests directory if it doesn't exist
mkdir -p tests

# Fix files with colons in names (if they exist)
if [ -f "test:test_core_function.py" ]; then
    mv "test:test_core_function.py" tests/test_core_functions.py
    echo "✓ Moved test:test_core_function.py → tests/test_core_functions.py"
fi

if [ -f "tests:confest.py" ]; then
    mv "tests:confest.py" tests/conftest.py
    echo "✓ Moved tests:confest.py → tests/conftest.py"
fi

if [ -f "tests:test_data_generator.py" ]; then
    mv "tests:test_data_generator.py" tests/test_data_generator.py
    echo "✓ Moved tests:test_data_generator.py → tests/test_data_generator.py"
fi

# Rename test.ini to pytest.ini
if [ -f "test.ini" ]; then
    mv test.ini pytest.ini
    echo "✓ Renamed test.ini → pytest.ini"
fi

# Move test_installation.py to tests/
if [ -f "test_installation.py" ]; then
    mv test_installation.py tests/
    echo "✓ Moved test_installation.py → tests/"
fi

# Create tests/__init__.py if it doesn't exist
touch tests/__init__.py
echo "✓ Created tests/__init__.py"

# Check if test_quick.py exists, if not create it
if [ ! -f "test_quick.py" ]; then
    echo "⚠️  test_quick.py not found. Creating it..."
    cat > test_quick.py << 'EOF'
#!/usr/bin/env python
"""Quick test to verify the analysis works with test data."""

import sys
import os
import time

def run_quick_test():
    """Run a quick test with the test data."""
    
    print("="*60)
    print("LipidProteinAnalyzer - Quick Test")
    print("="*60)
    
    try:
        import time_series_analysis as tsa
        print("✓ Module imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import: {e}")
        return False
    
    # Check test data
    if not os.path.exists('test_data/test_system.psf'):
        print("✗ Test data not found!")
        return False
    
    print("✓ Test data found")
    
    # Configure for test data
    tsa.TOPOLOGY_FILE = 'test_data/test_system.psf'
    tsa.TRAJECTORY_FILE = 'test_data/test_trajectory.xtc'
    tsa.OUTPUT_DIR = 'test_output'
    os.makedirs('test_output', exist_ok=True)
    
    try:
        # Quick test
        print("\nRunning quick analysis...")
        u = tsa.load_universe()
        leaflet0, leaflet1, L = tsa.identify_lipid_leaflets(u)
        lipid_sels = tsa.setup_lipid_selections(leaflet0, leaflet1)
        proteins = tsa.select_proteins(u, n_proteins=4)
        
        # Analyze just 5 frames
        time_series_data = tsa.analyze_time_series(
            u, proteins, lipid_sels, start=0, stop=5, step=1, interval=1
        )
        
        print(f"✓ Analysis completed: {len(time_series_data)} frames")
        print("\n✓ TEST PASSED!")
        return True
        
    except Exception as e:
        print(f"✗ TEST FAILED: {e}")
        return False

if __name__ == "__main__":
    success = run_quick_test()
    sys.exit(0 if success else 1)
EOF
    echo "✓ Created test_quick.py"
fi

echo ""
echo "File structure cleaned! Current structure:"
echo ""
echo "Root directory:"
ls -1 *.py *.md *.txt *.ini 2>/dev/null | head -10

echo ""
echo "tests/ directory:"
ls -1 tests/

echo ""
echo "test_data/ directory:"
ls -1 test_data/

echo ""
echo "✓ Cleanup complete!"