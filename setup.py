#!/usr/bin/env python3
"""
CIDD: Collaborative Intelligence for Structure-Based Drug Design Empowered by LLMs
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Read requirements
with open(os.path.join(this_directory, 'requirements.txt'), encoding='utf-8') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name='cidd',
    version='1.0.0',
    description='Collaborative Intelligence for Structure-Based Drug Design Empowered by LLMs',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Bowen Gao',
    author_email='bowen-gao@example.com',  # Please update with actual email
    url='https://github.com/bowen-gao/CIDD',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=requirements,
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
    keywords='drug design, molecular generation, LLM, structure-based drug design',
    project_urls={
        'Paper': 'https://openreview.net/pdf?id=7k7cubl1iL',
        'Source': 'https://github.com/bowen-gao/CIDD',
        'Tracker': 'https://github.com/bowen-gao/CIDD/issues',
    },
)