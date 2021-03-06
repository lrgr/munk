from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

import sys

if sys.version_info.major != 3:
        raise RuntimeError('MUNK requires Python 3')

setup(
    name='munk',
    version='0.1',

    description='MUNK: MUlti-species Network Kernel',
    url='https://github.com/lrgr/munk',
    author='Leiserson Research Group',

    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    keywords='bioinformatics',
    packages=['munk'],
    install_requires=[
        'numpy>=1.14.0,<2.0.0',
        'scipy>=0.19.1,<1.0.0',
        'scikit-learn>=0.19.1,<1.0.0',
        'networkx>=2.0,<3.0',
        ],

    project_urls={  # Optional
        'LRGR': 'https://lrgr.io',
    },
)
