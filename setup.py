# coding: utf-8
from __future__ import print_function, absolute_import
import sys
from os.path import join as pjoin

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

from presto import __version__, __author__, __license__

# from http://stackoverflow.com/questions/14399534/how-can-i-reference-requirements-txt-for-the-install-requires-kwarg-in-setuptool
install_reqs = parse_requirements("requirements.txt", session=False)
reqs = [str(ir.req) for ir in install_reqs]

_presto_scripts = ['AlignSets.py', 'AssemblePairs.py', 'BuildConsensus.py',
                   'ClusterSets.py', 'CollapseSeq.py', 'ConvertHeaders.py',
                   'EstimateError.py', 'FilterSeq.py',
                   'MaskPrimers.py', 'PairSeq.py', 'ParseHeaders.py',
                   'ParseLog.py', 'SplitSeq.py']

presto_scripts = [pjoin('presto', s) for s in _presto_scripts]

setup(
    name='presto',
    version=__version__, # can also use __version__
    author=__author__,
    author_email='jason.vanderheiden@yale.edu',
    description='A bioinformatics library for processing RepSeq data.',
    zip_safe=False,
    license=__license__,
    url='https://clip.med.yale.edu/presto',
    install_requires=reqs,
    packages=['presto'],
    scripts=presto_scripts,
    package_data={'': ['*.sh']},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)