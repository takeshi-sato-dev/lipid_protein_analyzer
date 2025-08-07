#!/usr/bin/env python
"""
Prepare repository for GitHub upload and JOSS submission.
This script helps ensure everything is ready.
"""

import os
import subprocess
from pathlib import Path
import json

def check_git_status():
    """Check if git is initialized and status."""
    try:
        result = subprocess.run(['git', 'status'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✅ Git repository already initialized")
            return True
        else:
            print("⚠️  Git not initialized")
            return False
    except FileNotFoundError:
        print("❌ Git is not installed!")
        return False

def create_gitignore():
    """Create a comprehensive .gitignore file."""
    gitignore_content = """# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
env/
venv/
.venv/
.env

# Testing
.coverage
.pytest_cache/
htmlcov/
*.cover
.hypothesis/
.tox/

# Jupyter Notebook
.ipynb_checkpoints

# Output directories
test_output/
test_data/
paper_figures/
paper_figures_test/
lipid_protein_timeseries/
results/
output/

# Large MD files
*.xtc
*.trr
*.dcd
*.pdb
*.gro
# Keep small test files
!test_data/*.psf
!test_data/*.xtc

# Temporary files
*.tmp
*.bak
*.backup
*~

# IDE
.vscode/
.idea/
*.swp
*.swo

# OS
.DS_Store
Thumbs.db
desktop.ini

# Documentation builds
docs/_build/
docs/_static/
docs/_templates/
"""
    
    gitignore_path = Path('.gitignore')
    gitignore_path.write_text(gitignore_content)
    print("✅ Created .gitignore")

def check_required_files():
    """Check all required files exist."""
    required_files = {
        'README.md': 'Main documentation',
        'LICENSE': 'MIT License',
        'requirements.txt': 'Python dependencies',
        'setup.py': 'Installation script',
        'paper.md': 'JOSS paper',
        'paper.bib': 'Bibliography',
        'time_series_analysis.py': 'Main code',
        'test_quick.py': 'Quick test',
        'run_analysis.py': 'Analysis runner'
    }
    
    print("\n📋 Checking required files:")
    all_present = True
    
    for file, description in required_files.items():
        if Path(file).exists():
            print(f"  ✅ {file:30} - {description}")
        else:
            print(f"  ❌ {file:30} - MISSING!")
            all_present = False
    
    return all_present

def update_readme_template():
    """Add GitHub-specific badges to README."""
    readme_path = Path('README.md')
    if not readme_path.exists():
        print("❌ README.md not found")
        return
    
    readme_content = readme_path.read_text()
    
    # Check if badges already exist
    if '![Tests]' in readme_content or '[![DOI]' in readme_content:
        print("✅ README already has badges")
        return
    
    # Add badges at the beginning
    badges = """# LipidProteinAnalyzer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
![Tests](https://github.com/YOUR_USERNAME/lipid_protein_analyzer/workflows/tests/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

"""
    
    if readme_content.startswith('# LipidProteinAnalyzer'):
        readme_content = readme_content.replace('# LipidProteinAnalyzer\n', badges)
    else:
        readme_content = badges + readme_content
    
    # Save backup
    backup_path = Path('README.md.backup')
    readme_path.rename(backup_path)
    
    # Write updated version
    readme_path.write_text(readme_content)
    print("✅ Updated README with badge placeholders")
    print("   ⚠️  Remember to update YOUR_USERNAME and ZENODO DOI after setup")

def create_citation_file():
    """Create CITATION.cff file for GitHub."""
    citation_content = """cff-version: 1.2.0
message: "If you use this software, please cite it as below."
authors:
  - family-names: "Your Last Name"
    given-names: "Your First Name"
    orcid: "https://orcid.org/0009-0006-9156-8655"
title: "LipidProteinAnalyzer: Time-series analysis of lipid-protein interactions"
version: 1.0.0
date-released: 2025-01-07
url: "https://github.com/YOUR_USERNAME/lipid_protein_analyzer"
doi: 10.5281/zenodo.XXXXXXX
license: MIT
keywords:
  - molecular dynamics
  - lipid-protein interactions
  - membrane proteins
  - time series analysis
  - parallel processing
"""
    
    citation_path = Path('CITATION.cff')
    citation_path.write_text(citation_content)
    print("✅ Created CITATION.cff")
    print("   ⚠️  Remember to update author information and DOI")

def setup_github_actions():
    """Create GitHub Actions workflow for CI."""
    workflows_dir = Path('.github/workflows')
    workflows_dir.mkdir(parents=True, exist_ok=True)
    
    # Simple test workflow
    workflow_content = """name: Tests

on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ['3.8', '3.10']
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install pytest
    
    - name: Run tests
      run: |
        pytest tests/ -v
"""
    
    workflow_path = workflows_dir / 'tests.yml'
    workflow_path.write_text(workflow_content)
    print("✅ Created GitHub Actions workflow")

def initialize_git():
    """Initialize git repository."""
    if not check_git_status():
        print("\n🔧 Initializing Git repository...")
        subprocess.run(['git', 'init'])
        subprocess.run(['git', 'branch', '-M', 'main'])
        print("✅ Git repository initialized")

def main():
    print("="*60)
    print("🚀 Preparing for GitHub Upload and JOSS Submission")
    print("="*60)
    
    # Check required files
    if not check_required_files():
        print("\n⚠️  Some required files are missing!")
        print("Please ensure all files are present before proceeding.")
        return 1
    
    # Create/update files
    print("\n📝 Creating/updating files...")
    create_gitignore()
    update_readme_template()
    create_citation_file()
    setup_github_actions()
    
    # Initialize git
    initialize_git()
    
    print("\n" + "="*60)
    print("✅ Repository prepared for GitHub!")
    print("="*60)
    
    print("\n📋 Next steps:")
    print("1. Create repository on GitHub:")
    print("   https://github.com/new")
    print("   Name: lipid_protein_analyzer")
    print("   Visibility: Public")
    print("")
    print("2. Add all files to git:")
    print("   git add .")
    print('   git commit -m "Initial commit for JOSS submission"')
    print("")
    print("3. Add GitHub remote (replace YOUR_USERNAME):")
    print("   git remote add origin https://github.com/YOUR_USERNAME/lipid_protein_analyzer.git")
    print("   git push -u origin main")
    print("")
    print("4. Set up Zenodo:")
    print("   - Go to https://zenodo.org/account/settings/github/")
    print("   - Enable your repository")
    print("   - Create a release on GitHub")
    print("   - Get DOI from Zenodo")
    print("")
    print("5. Update files with actual values:")
    print("   - Replace YOUR_USERNAME in README.md")
    print("   - Replace XXXXXXX with actual Zenodo DOI")
    print("   - Update author name in CITATION.cff")
    print("   - Update author name in paper.md")
    print("")
    print("6. Submit to JOSS:")
    print("   https://joss.theoj.org/papers/new")
    
    return 0

if __name__ == "__main__":
    exit(main())