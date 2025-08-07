"""Setup script for LipidProteinAnalyzer using original time_series_analysis.py"""

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="lipid-protein-analyzer",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Time-series analysis of lipid-protein interactions in MD simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/lipid_protein_analyzer",
    py_modules=["time_series_analysis"],  # あなたの元のファイル名
    python_requires=">=3.8",
    install_requires=[
        "MDAnalysis>=2.0.0",
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "tqdm>=4.60.0",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    # コマンドラインツールとして登録（オプション）
    entry_points={
        'console_scripts': [
            'lpa-analyze=time_series_analysis:main',
        ],
    },
)