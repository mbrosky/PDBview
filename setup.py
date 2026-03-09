"""Setup script for PDBView"""

from setuptools import setup, find_packages

setup(
    name="pdbview",
    version="0.2.0",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "textual>=0.50.0",
        "biopython>=1.83",
        "numpy>=1.24.0",
        "requests>=2.31.0",
    ],
    extras_require={
        "rdkit": ["rdkit>=2023.9.1"],
        "full": ["rdkit>=2023.9.1", "prolif>=2.0.0"],
    },
    entry_points={
        "console_scripts": [
            "pdbview=pdbview.cli:main",
        ],
    },
    author="Briki Mondher",
    author_email="mondher.briki.isdd@gmail.com",
    description="Terminal protein and ligand structure viewer",
    keywords="protein ligand pdb mol2 sdf terminal tui bioinformatics cheminformatics",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)
