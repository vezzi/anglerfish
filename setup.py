#!/usr/bin/env python
"""
Anglerfish is a tool designed to demultiplex Illumina libraries sequenced on Oxford Nanopore flowcells. The primary purpose for this would be to do QC, i.e. to check pool balancing, assess contamination, library insert sizes and so on.

Install with pip:

    pip install bio-anglerfish

Or from bioconda:

    conda install -c bioconda anglerfish
"""

from pathlib import Path

from setuptools import find_packages, setup

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

version = "0.7.0"


def _load_requirements(path: str = "requirements.txt") -> list[str]:
    reqs: list[str] = []
    with open(path) as f:
        for line in f:
            dep = line.strip()
            if dep and not dep.startswith("#"):
                reqs.append(dep)
    return reqs

setup(
    name="bio-anglerfish",
    version=version,
    description="Anglerfish, a tool to demultiplex Illumina libraries from ONT data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Remi-Andre Olsen",
    author_email="remi-andre.olsen@scilifelab.se",
    url="https://github.com/NationalGenomicsInfrastructure/anglerfish",
    license="MIT",
    python_requires=">=3.12",
    packages=find_packages(),
    package_data={"": ["config/adaptors.yaml"]},
    install_requires=_load_requirements(),
    extras_require={
        "dev": [
            "pytest>=8.4",
            "ruff>=0.15",
            "mypy>=1.19",
            "editorconfig-checker>=3.6",
            "pre-commit>=4.3",
        ]
    },
    entry_points={
        "console_scripts": [
            "anglerfish=anglerfish.cli:app",
        ],
    },
    zip_safe=False,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
