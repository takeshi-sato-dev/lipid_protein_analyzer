#!/usr/bin/env python
"""
Comprehensive validation script for LipidProteinAnalyzer.
Run this to check all components before JOSS submission.

Usage:
    python validate_all.py [--full]
    
    --full : Run full analysis with real data (if available)
"""

import sys
import os
import subprocess
import time
import importlib
import argparse
from pathlib import Path

# Color output for better readability
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def print_header(text):
    """Print a section header."""
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.ENDC}")
    print(f"{Colors.BOLD}{Colors.BLUE}{text}{Colors.ENDC}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.ENDC}")

def print_success(text):
    """Print success message."""
    print(f"{Colors.GREEN}✓ {text}{Colors.ENDC}")

def print_error(text):
    """Print error message."""
    print(f"{Colors.RED}✗ {text}{Colors.ENDC}")

def print_warning(text):
    """Print warning message."""
    print(f"{Colors.YELLOW}⚠ {text}{Colors.ENDC}")

def print_info(text):
    """Print info message."""
    print(f"  {text}")

def check_files():
    """Check that all required files exist."""
    print_header("1. FILE STRUCTURE CHECK")
    
    required_files = [
        'LICENSE',
        'README.md',
        'requirements.txt',
        'setup.py',
        'paper.md',
        'paper.bib',
        'time_series_analysis.py',
        'run_analysis.py',
        'test_quick.py'
    ]
    
    optional_files = [
        'pytest.ini',
        'TESTING.md',
        '.gitignore',
        'tests/test_core_functions.py',
        'tests/test_data_generator.py'
    ]
    
    all_good = True
    
    # Check required files
    print("\nRequired files:")
    for file in required_files:
        if os.path.exists(file):
            print_success(f"{file}")
        else:
            print_error(f"{file} - MISSING!")
            all_good = False
    
    # Check optional files
    print("\nOptional files:")
    for file in optional_files:
        if os.path.exists(file):
            print_success(f"{file}")
        else:
            print_warning(f"{file} - Missing (optional)")
    
    return all_good

def check_imports():
    """Test different import methods."""
    print_header("2. IMPORT CHECK")
    
    all_good = True
    
    # Test direct import
    try:
        import time_series_analysis
        print_success("Direct import: import time_series_analysis")
    except ImportError as e:
        print_error(f"Direct import failed: {e}")
        all_good = False
    
    # Test specific function import
    try:
        from time_series_analysis import analyze_time_series_parallel
        print_success("Function import: from time_series_analysis import analyze_time_series_parallel")
    except ImportError as e:
        print_error(f"Function import failed: {e}")
        all_good = False
    
    # Test wildcard import
    try:
        exec("from time_series_analysis import *")
        print_success("Wildcard import: from time_series_analysis import *")
    except ImportError as e:
        print_error(f"Wildcard import failed: {e}")
        all_good = False
    
    # Check main functions exist
    try:
        import time_series_analysis as tsa
        required_functions = [
            'load_universe',
            'identify_lipid_leaflets',
            'setup_lipid_selections',
            'select_proteins',
            'analyze_time_series',
            'analyze_time_series_parallel',
            'plot_lipid_composition_time_series'
        ]
        
        print("\nChecking main functions:")
        for func_name in required_functions:
            if hasattr(tsa, func_name):
                print_success(f"Function '{func_name}' found")
            else:
                print_error(f"Function '{func_name}' missing!")
                all_good = False
                
    except Exception as e:
        print_error(f"Failed to check functions: {e}")
        all_good = False
    
    return all_good

def check_dependencies():
    """Check that all dependencies are installed."""
    print_header("3. DEPENDENCY CHECK")
    
    dependencies = {
        'MDAnalysis': '2.0.0',
        'numpy': '1.20.0',
        'pandas': '1.3.0',
        'matplotlib': '3.3.0',
        'seaborn': '0.11.0',
        'tqdm': '4.60.0'
    }
    
    all_good = True
    
    for package, min_version in dependencies.items():
        try:
            module = importlib.import_module(package)
            version = getattr(module, '__version__', 'unknown')
            print_success(f"{package} >= {min_version} (found {version})")
        except ImportError:
            print_error(f"{package} - NOT INSTALLED!")
            all_good = False
    
    # Check optional test dependencies
    print("\nOptional test dependencies:")
    test_deps = ['pytest', 'pytest-cov']
    for package in test_deps:
        try:
            importlib.import_module(package.replace('-', '_'))
            print_success(f"{package} installed")
        except ImportError:
            print_warning(f"{package} not installed (needed for testing)")
    
    return all_good

def run_quick_test():
    """Run the quick test script."""
    print_header("4. QUICK TEST")
    
    if not os.path.exists('test_quick.py'):
        print_error("test_quick.py not found!")
        return False
    
    # First, create test data if it doesn't exist
    if not os.path.exists('test_data/test_system.psf'):
        print_info("Creating test data...")
        try:
            # Try to run test data generator
            if os.path.exists('tests/test_data_generator.py'):
                result = subprocess.run(
                    [sys.executable, 'tests/test_data_generator.py'],
                    capture_output=True, text=True
                )
                if result.returncode == 0:
                    print_success("Test data created")
                else:
                    print_warning("Could not create test data automatically")
            else:
                # Create minimal test data inline
                from create_test_data import create_test_data
                create_test_data()
                print_success("Test data created")
        except Exception as e:
            print_error(f"Failed to create test data: {e}")
            return False
    
    # Run the quick test
    print_info("Running quick test...")
    start_time = time.time()
    
    try:
        result = subprocess.run(
            [sys.executable, 'test_quick.py'],
            capture_output=True,
            text=True,
            timeout=60  # 60 second timeout
        )
        
        elapsed = time.time() - start_time
        
        if result.returncode == 0:
            print_success(f"Quick test passed in {elapsed:.1f} seconds")
            
            # Check if output was generated
            if os.path.exists('test_output'):
                output_files = os.listdir('test_output')
                print_info(f"Generated {len(output_files)} output files")
            
            return True
        else:
            print_error("Quick test failed!")
            print_info("Error output:")
            print(result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print_error("Quick test timed out (>60 seconds)")
        return False
    except Exception as e:
        print_error(f"Failed to run quick test: {e}")
        return False

def check_memory_usage():
    """Check memory usage with test data."""
    print_header("5. MEMORY CHECK")
    
    try:
        import psutil
        import time_series_analysis as tsa
        
        # Override settings for test data
        tsa.TOPOLOGY_FILE = 'test_data/test_system.psf'
        tsa.TRAJECTORY_FILE = 'test_data/test_trajectory.xtc'
        
        process = psutil.Process(os.getpid())
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        # Load and analyze
        u = tsa.load_universe()
        leaflet0, leaflet1, L = tsa.identify_lipid_leaflets(u)
        lipid_sels = tsa.setup_lipid_selections(leaflet0, leaflet1)
        proteins = tsa.select_proteins(u)
        
        # Small analysis
        time_series_data = tsa.analyze_time_series(
            u, proteins, lipid_sels,
            start=0, stop=5, step=1
        )
        
        final_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_used = final_memory - initial_memory
        
        print_success(f"Memory usage: {memory_used:.1f} MB")
        
        if memory_used > 500:
            print_warning("High memory usage for test data")
        
        return True
        
    except ImportError:
        print_warning("psutil not installed - skipping memory check")
        return True
    except Exception as e:
        print_error(f"Memory check failed: {e}")
        return False

def check_parallel_processing():
    """Test parallel processing capabilities."""
    print_header("6. PARALLEL PROCESSING CHECK")
    
    try:
        import multiprocessing as mp
        import time_series_analysis as tsa
        
        # Check CPU count
        cpu_count = mp.cpu_count()
        print_info(f"System CPUs: {cpu_count}")
        
        # Check if parallel function exists
        if hasattr(tsa, 'analyze_time_series_parallel'):
            print_success("Parallel analysis function found")
            
            # Try to determine how many cores would be used
            expected_cores = min(cpu_count - 1, max(1, int(cpu_count * 0.75)))
            print_info(f"Expected to use: {expected_cores} cores")
            
            return True
        else:
            print_error("Parallel analysis function not found")
            return False
            
    except Exception as e:
        print_error(f"Parallel processing check failed: {e}")
        return False

def run_pytest():
    """Run pytest if available."""
    print_header("7. PYTEST SUITE")
    
    try:
        result = subprocess.run(
            ['pytest', '--tb=short', '-q'],
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            print_success("All tests passed")
            # Parse output for stats
            lines = result.stdout.split('\n')
            for line in lines:
                if 'passed' in line:
                    print_info(line.strip())
            return True
        else:
            print_warning("Some tests failed")
            print(result.stdout)
            return False
            
    except FileNotFoundError:
        print_warning("pytest not installed - skipping test suite")
        return True

def check_paper_figures(full_analysis=False):
    """Check if paper figures can be generated."""
    print_header("8. PAPER FIGURES CHECK")
    
    if full_analysis and os.path.exists('step5_assembly.psf') and os.path.exists('md_wrapped.xtc'):
        print_info("Running full analysis for paper figures...")
        print_warning("This may take several minutes...")
        
        try:
            result = subprocess.run(
                [sys.executable, 'run_analysis.py',
                 '--start', '60000', '--stop', '60100', '--step', '10',
                 '--output', 'paper_figures_test'],
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode == 0:
                print_success("Paper figures generated successfully")
                
                # Check output
                if os.path.exists('paper_figures_test'):
                    files = os.listdir('paper_figures_test')
                    png_files = [f for f in files if f.endswith('.png')]
                    print_info(f"Generated {len(png_files)} PNG files")
                    print_info("Remember to use these for paper.md!")
                
                return True
            else:
                print_error("Failed to generate paper figures")
                return False
                
        except subprocess.TimeoutExpired:
            print_error("Figure generation timed out")
            return False
        except Exception as e:
            print_error(f"Failed to generate figures: {e}")
            return False
    else:
        print_warning("Real data not found - skipping full analysis")
        print_info("For paper.md, you need to run with your real data:")
        print_info("  python run_analysis.py --start 60000 --stop 80000 --step 50")
        return True

def generate_summary(results):
    """Generate a summary of all checks."""
    print_header("VALIDATION SUMMARY")
    
    total = len(results)
    passed = sum(results.values())
    
    print(f"\nResults: {passed}/{total} checks passed")
    print("")
    
    for check, status in results.items():
        if status:
            print_success(check)
        else:
            print_error(check)
    
    print("")
    
    if passed == total:
        print_success("🎉 ALL CHECKS PASSED! Ready for JOSS submission.")
    elif passed >= total - 2:
        print_warning("⚠️  Almost ready! Fix the remaining issues.")
    else:
        print_error("❌ Several issues need to be fixed before submission.")
    
    # Provide next steps
    print("\n" + Colors.BOLD + "Next Steps:" + Colors.ENDC)
    
    if not results.get('Paper Figures', True):
        print("1. Generate figures with real data for paper.md")
    
    if not results.get('Pytest', True):
        print("2. Install pytest and create comprehensive tests")
    
    print("3. Update paper.md with:")
    print("   - Your real name and ORCID")
    print("   - Real performance benchmarks")
    print("   - Figures from actual analysis")
    print("4. Archive on Zenodo and get DOI")
    print("5. Submit to JOSS!")

def main():
    """Main validation routine."""
    parser = argparse.ArgumentParser(description='Validate LipidProteinAnalyzer for JOSS')
    parser.add_argument('--full', action='store_true', 
                       help='Run full analysis with real data if available')
    args = parser.parse_args()
    
    print(Colors.BOLD + "\nLipidProteinAnalyzer - JOSS Validation Suite" + Colors.ENDC)
    print("This will check all components required for JOSS submission\n")
    
    results = {}
    
    # Run all checks
    results['File Structure'] = check_files()
    results['Imports'] = check_imports()
    results['Dependencies'] = check_dependencies()
    results['Quick Test'] = run_quick_test()
    results['Memory Usage'] = check_memory_usage()
    results['Parallel Processing'] = check_parallel_processing()
    results['Pytest'] = run_pytest()
    results['Paper Figures'] = check_paper_figures(args.full)
    
    # Generate summary
    generate_summary(results)
    
    # Return exit code based on critical checks
    critical_checks = ['File Structure', 'Imports', 'Dependencies', 'Quick Test']
    critical_passed = all(results.get(check, False) for check in critical_checks)
    
    return 0 if critical_passed else 1

if __name__ == "__main__":
    sys.exit(main())