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
    version='1.3.1',
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
    install_requires=[
        'numpy>=1.16.4',
        'pandas>=0.24.2',
        'tables>=3.5.2',
        'scipy>=1.3.0',
        'matplotlib>=3.1.0',
        'scikit-learn>=0.21.2',
        'plotly>=4.1.1',
        'mygene>=3.1.0',
        'pynndescent>=0.3.3',
        'networkx>=2.3',
        'python-louvain>=0.13',
        'adjustText>=0.7.3',
        'umap-learn>=0.3.10',
        'phate>=1.0.3',
        'fitsne>=1.0.1; platform_system=="Linux" or platform_system=="Darwin"'],
    zip_safe=False
)
