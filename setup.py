import os
from setuptools import Extension, setup, find_packages
from Cython.Build import cythonize

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    ext_modules=cythonize("src/ChromBridGE/ChromBridGE_aln.pyx"),
    name = "ChromBridGE",
    version = "0.0.5",
    author = "Kendell Clement",
    author_email = "kclement@mgh.harvard.edu",
    description = "Detection of CRISPR-Cas-associated translocations from next-generation sequencing data",
    long_description=long_description,
    keywords = "NGS Translocations",
    package_dir={"": "src"},
    packages=find_packages(where='src'),
    py_modules=['ChromBridGE','ChromBridGE.ChromBridGE_tx'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
    ],
    entry_points={
        'console_scripts': [
            'ChromBridGE = ChromBridGE.ChromBridGE:main',
        ]
    },
    install_requires=[
        'numpy',
        'Cython',
    ],
)
