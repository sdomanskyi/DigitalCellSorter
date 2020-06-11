"""
setup.py for DigitalCellSorter package
"""
from setuptools import setup, find_packages
from codecs import open
from os import path
here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description=f.read()

setup(
    name='DigitalCellSorter',
    packages=find_packages(),
    version='1.3.4.10',
    description='Toolkit for analysis and identification of cell types from heterogeneous single cell RNA-seq data',
    long_description_content_type="text/markdown",
    long_description=long_description,
    include_package_data=True,
    author='S. Domanskyi , A. Szedlak, N. T Hawkins, J. Wang, G. Paternostro, C. Piermarocchi',
    author_email='s.domanskyi@gmail.com',
    license='MIT License',
    url='https://github.com/sdomanskyi/DigitalCellSorter',
    download_url='https://github.com/sdomanskyi/DigitalCellSorter/archive/1.3.0.tar.gz',
    keywords=['single cell RNA sequencing', 'cell type identification','biomarkers'],
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: Unix',
        'Topic :: Education',
        'Topic :: Utilities',
        ],
    python_requires='>=3',
    install_requires=[
        'numpy>=1.16.4',
        'pandas>=0.24.2',
        'patsy>=0.5.1',
        'xlrd>=1.2.0',
        'openpyxl>=3.0.3',
        'tables>=3.4.2',
        'scipy>=1.3.0',
        'matplotlib>=3.1.0',
        'scikit-learn>=0.21.2',
        'mygene>=3.1.0',
        'plotly>=4.1.1',
        'adjustText>=0.7.3'],
    zip_safe=False
)
