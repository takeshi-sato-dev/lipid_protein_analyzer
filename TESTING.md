# Testing Guide

## Quick Start

```bash
# Install test dependencies
pip install pytest pytest-cov

# Generate test data
python tests/test_data_generator.py

# Run all tests
pytest

# Run with coverage
pytest --cov=time_series_analysis --cov-report=html
```

## Test Structure

```
tests/
├── test_core_functions.py    # Core functionality tests (~250 lines)
├── test_data_generator.py    # Generate minimal test data (~100 lines)
├── test_data/                # Generated test files (auto-created)
│   ├── test_system.psf      # Minimal topology (1 KB)
│   └── test_trajectory.xtc  # 10 frames (10 KB)
└── conftest.py              # Pytest fixtures (optional)
```

## Test Categories

- **Unit Tests**: Test individual functions with mock data
- **Integration Tests**: Test complete workflows with small test data
- **Visualization Tests**: Verify plot generation without large datasets

## Running Specific Tests

```bash
# Run only fast tests
pytest -m "not slow"

# Run specific test file
pytest tests/test_core_functions.py

# Run with verbose output
pytest -v

# Run specific test class
pytest tests/test_core_functions.py::TestFrameProcessing
```

## Continuous Integration

Tests automatically run on GitHub Actions for:
- Python 3.8, 3.9, 3.10, 3.11
- Ubuntu and macOS
- Coverage reports uploaded to Codecov

## JOSS Requirements

✅ Automated tests that cover core functionality  
✅ Tests run with small data (<1 MB total)  
✅ CI/CD integration with GitHub Actions  
✅ Clear documentation for running tests