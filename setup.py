#!/usr/bin/env python
"""
Presto setup
"""
# Future
from __future__ import division, absolute_import, print_function

# Imports
import os, re, sys

# Check setup requirements
if sys.version_info < (2,7,5):
    print('At least Python 2.7.5 is required.\n', file=sys.stderr)
    exit(1)

try:
    from setuptools import setup
except ImportError:
    print("Please install setuptools before installing presto.", file=sys.stderr)
    exit(1)

try:
    from pip.req import parse_requirements
except ImportError:
    print("Please install pip before installing presto.", file=sys.stderr)
    exit(1)


# Get version
# Modified from requests setup
version = ''
with open('presto/_version.py', 'r') as handle:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        handle.read(), re.MULTILINE).group(1)
if not version:
    raise RuntimeError('Cannot find version information')


# Parse requirements
requirements = parse_requirements("requirements.txt", session=False)
install_requires = [str(r.req) for r in requirements]

# Define installation path for commandline tools
scripts = ['AlignSets.py',
           'AssemblePairs.py',
           'BuildConsensus.py',
           'ClusterSets.py',
           'CollapseSeq.py',
           'ConvertHeaders.py',
           'EstimateError.py',
           'FilterSeq.py',
           'MaskPrimers.py',
           'PairSeq.py',
           'ParseHeaders.py',
           'ParseLog.py',
           'SplitSeq.py']
install_scripts = [os.path.join('bin', s) for s in scripts]

# Load long package description
with open('README.md', 'r') as f:
    long_description = ''.join([x for x in f])

# Setup
setup(name='presto',
      version=version,
      author='Jason Vander Heiden',
      author_email='jason.vanderheiden@yale.edu',
      description='A bioinformatics toolkit for processing high-throughput lymphocyte receptor sequencing data.',
      long_description=long_description,
      zip_safe=False,
      license='Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported',
      url='https://clip.med.yale.edu/presto',
      keywords='bioinformatics immunoglobulin lymphocyte sequencing',
      install_requires=install_requires,
      packages=['presto'],
      scripts=install_scripts,
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'])