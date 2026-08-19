"""
Setup script for leukoquant package.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# Read requirements
requirements = (this_directory / "requirements.txt").read_text().splitlines()

setup(
    name="leukoquant",
    version="1.0.0",
    author="Stylianos Charalampous",
    author_email="stylianosc@users.noreply.github.com",
    description="Lesion-informed white matter damage metrics toolkit for cerebral small vessel disease",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/stylianosc/leukoquant",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Image Processing",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
        ]
    },
    entry_points={
        "console_scripts": [
            "leukoquant=leukoquant.cli.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "leukoquant": [
            "mappings/*",
            "external/**/*",
        ]
    },
    zip_safe=False,
)
