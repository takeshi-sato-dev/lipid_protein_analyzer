name: Tests

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main ]
  workflow_dispatch:  # Allow manual triggering

jobs:
  test:
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest]
        python-version: ['3.8', '3.9', '3.10', '3.11']
        exclude:
          # Skip some combinations to save CI time
          - os: macos-latest
            python-version: '3.8'
          - os: macos-latest
            python-version: '3.9'
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}
    
    - name: Cache pip packages
      uses: actions/cache@v3
      with:
        path: ~/.cache/pip
        key: ${{ runner.os }}-pip-${{ hashFiles('requirements.txt') }}
        restore-keys: |
          ${{ runner.os }}-pip-
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install pytest pytest-cov pytest-timeout
        pip install -e .
    
    - name: Generate test data
      run: |
        python tests/test_data_generator.py
    
    - name: Run quick test
      run: |
        python test_quick.py
      timeout-minutes: 2
    
    - name: Run pytest with coverage
      run: |
        pytest tests/ \
          --cov=time_series_analysis \
          --cov-report=xml \
          --cov-report=term-missing \
          --timeout=60 \
          -v
      timeout-minutes: 10
    
    - name: Upload coverage to Codecov
      if: matrix.python-version == '3.10' && matrix.os == 'ubuntu-latest'
      uses: codecov/codecov-action@v3
      with:
        file: ./coverage.xml
        flags: unittests
        name: codecov-umbrella
        fail_ci_if_error: false

  validation:
    runs-on: ubuntu-latest
    needs: test
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.10'
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install -e .
    
    - name: Create test data
      run: |
        python create_test_data.py
    
    - name: Run validation script
      run: |
        python validate_all.py
      continue-on-error: true  # Don't fail CI if validation has warnings
    
    - name: Check file structure
      run: |
        # Check all required files exist
        for file in LICENSE README.md requirements.txt setup.py paper.md paper.bib; do
          if [ ! -f "$file" ]; then
            echo "Missing required file: $file"
            exit 1
          fi
        done
        echo "All required files present ✓"

  performance:
    runs-on: ubuntu-latest
    if: github.event_name == 'push' && github.ref == 'refs/heads/main'
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.10'
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install memory_profiler
        pip install -e .
    
    - name: Generate test data
      run: |
        python tests/test_data_generator.py
    
    - name: Profile memory usage
      run: |
        python -m memory_profiler test_quick.py > memory_profile.txt 2>&1 || true
        
        # Extract peak memory if available
        if grep -q "MiB" memory_profile.txt; then
          echo "Memory profiling results:"
          grep "MiB" memory_profile.txt | head -5
        fi
      continue-on-error: true
    
    - name: Test parallel performance
      run: |
        # Create a simple performance test script
        cat << 'EOF' > test_performance.py
        import time
        import os
        os.environ['TESTING'] = '1'  # Skip any GUI elements
        
        from time_series_analysis import *
        
        # Load test data
        u = load_universe()
        leaflet0, leaflet1, L = identify_lipid_leaflets(u)
        lipid_sels = setup_lipid_selections(leaflet0, leaflet1)
        proteins = select_proteins(u)
        
        # Time serial execution
        start = time.time()
        serial_result = analyze_time_series(
            u, proteins, lipid_sels, 
            start=0, stop=5, step=1
        )
        serial_time = time.time() - start
        
        # Time parallel execution
        start = time.time()
        parallel_result = analyze_time_series_parallel(
            u, proteins, lipid_sels,
            start=0, stop=5, step=1
        )
        parallel_time = time.time() - start
        
        print(f"Serial time: {serial_time:.2f}s")
        print(f"Parallel time: {parallel_time:.2f}s")
        
        if parallel_time < serial_time:
            print(f"Speedup: {serial_time/parallel_time:.2f}x ✓")
        else:
            print("No speedup observed (expected for small test data)")
        EOF
        
        python test_performance.py
      continue-on-error: true

  docs:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Check documentation
      run: |
        # Check README has essential sections
        for section in "Installation" "Usage" "Features" "Citation"; do
          if ! grep -q "$section" README.md; then
            echo "README.md missing section: $section"
            exit 1
          fi
        done
        echo "README.md structure ✓"
        
        # Check paper.md is valid
        if ! grep -q "^title:" paper.md; then
          echo "paper.md missing title"
          exit 1
        fi
        if ! grep -q "^authors:" paper.md; then
          echo "paper.md missing authors"
          exit 1
        fi
        echo "paper.md structure ✓"