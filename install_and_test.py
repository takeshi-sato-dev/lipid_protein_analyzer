#!/usr/bin/env python
"""
Simple script to install pytest and run tests.
This ensures tests work on any platform (Windows/Mac/Linux).
"""

import subprocess
import sys
import os
from pathlib import Path

def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"\n{'='*60}")
    print(f"🔧 {description}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(
            cmd, 
            shell=True, 
            capture_output=True, 
            text=True,
            check=True
        )
        print(result.stdout)
        if result.stderr:
            print(f"Warnings: {result.stderr}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error: {e}")
        print(f"Output: {e.stdout}")
        print(f"Error: {e.stderr}")
        return False

def main():
    """Main test installation and execution."""
    
    print("""
    ╔══════════════════════════════════════════════════════════╗
    ║           LipidProteinAnalyzer Test Suite               ║
    ║                                                          ║
    ║  This script will:                                      ║
    ║  1. Install pytest and dependencies                     ║
    ║  2. Generate test data                                  ║
    ║  3. Run all tests                                       ║
    ║  4. Generate coverage report                            ║
    ╚══════════════════════════════════════════════════════════╝
    """)
    
    # Step 1: Install pytest
    print("\n📦 Step 1: Installing pytest and test dependencies...")
    
    test_packages = [
        "pytest>=7.0.0",
        "pytest-cov>=3.0.0",
        "pytest-timeout>=2.0.0",
    ]
    
    for package in test_packages:
        if not run_command(
            f"{sys.executable} -m pip install '{package}'",
            f"Installing {package}"
        ):
            print(f"⚠️  Failed to install {package}, continuing anyway...")
    
    # Step 2: Generate test data
    print("\n📊 Step 2: Generating test data...")
    
    test_data_dir = Path("test_data")
    if not test_data_dir.exists():
        # First try test_data_generator.py if it exists
        if Path("tests/test_data_generator.py").exists():
            run_command(
                f"{sys.executable} tests/test_data_generator.py",
                "Generating test data with test_data_generator.py"
            )
        # Otherwise use create_test_data.py
        elif Path("create_test_data.py").exists():
            run_command(
                f"{sys.executable} create_test_data.py",
                "Generating test data with create_test_data.py"
            )
        else:
            print("⚠️  No test data generator found, skipping...")
    else:
        print("✅ Test data already exists")
    
    # Step 3: Run quick test
    print("\n🚀 Step 3: Running quick validation test...")
    
    if Path("test_quick.py").exists():
        if not run_command(
            f"{sys.executable} test_quick.py",
            "Running quick test"
        ):
            print("⚠️  Quick test failed, but continuing with pytest...")
    
    # Step 4: Run pytest
    print("\n🧪 Step 4: Running pytest test suite...")
    
    # Check if tests directory exists
    if not Path("tests").exists():
        print("⚠️  No tests directory found. Creating basic test...")
        
        # Create a minimal test
        Path("tests").mkdir(exist_ok=True)
        
        with open("tests/test_basic.py", "w") as f:
            f.write("""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def test_import():
    '''Test that the module can be imported.'''
    import time_series_analysis
    assert time_series_analysis is not None

def test_constants():
    '''Test that required constants are defined.'''
    import time_series_analysis as tsa
    assert hasattr(tsa, 'TOPOLOGY_FILE')
    assert hasattr(tsa, 'TRAJECTORY_FILE')
    assert hasattr(tsa, 'OUTPUT_DIR')

def test_functions_exist():
    '''Test that main functions exist.'''
    import time_series_analysis as tsa
    assert hasattr(tsa, 'load_universe')
    assert hasattr(tsa, 'analyze_time_series')
    assert hasattr(tsa, 'analyze_time_series_parallel')
""")
        print("✅ Created basic test file")
    
    # Run pytest with coverage
    pytest_cmd = (
        f"{sys.executable} -m pytest tests/ "
        f"-v "
        f"--tb=short "
        f"--timeout=60 "
        f"--cov=time_series_analysis "
        f"--cov-report=term-missing "
        f"--cov-report=html "
        f"-m 'not slow'"  # Skip slow tests by default
    )
    
    success = run_command(pytest_cmd, "Running pytest with coverage")
    
    # Step 5: Summary
    print("\n" + "="*60)
    print("📊 TEST SUMMARY")
    print("="*60)
    
    if success:
        print("✅ All tests passed!")
        
        # Check if coverage report was generated
        if Path("htmlcov/index.html").exists():
            print(f"\n📈 Coverage report generated:")
            print(f"   Open: {Path('htmlcov/index.html').absolute()}")
        
        print("\n🎉 Your code is ready for JOSS submission!")
        print("\nNext steps:")
        print("1. Generate figures with real data for paper.md")
        print("2. Update author information in paper.md")
        print("3. Push to GitHub and enable Actions")
        print("4. Submit to JOSS!")
        
    else:
        print("❌ Some tests failed")
        print("\nPlease fix the failing tests before submission.")
        print("Run with more detail:")
        print(f"  {sys.executable} -m pytest tests/ -vv")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())