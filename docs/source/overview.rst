**Overview**
============

Identification of hematological cell types from heterogeneous single cell RNA-seq data.

Polled Digital Cell Sorter (p-DCS): Automatic identification of hematological cell types from single cell RNA-sequencing clusters.
Sergii Domanskyi, Anthony Szedlak, Nathaniel T Hawkins, Jiayin Wang, Giovanni Paternostro & Carlo Piermarocchi, *BMC Bioinformatics* volume 20, Article number: 369 (**2019**) https://doi.org/10.1186/s12859-019-2951-x


**Getting Started**
===================

These instructions will get you a copy of the project up and running on your machine for data analysis, development or testing purposes.

**Installation**
----------------

The code runs in Python >= 3.7 environment. 

It is highly recommended to install Anaconda.
Installers are available at https://www.anaconda.com/distribution/

It uses packages ``numpy``, ``pandas``, ``matplotlib``, ``scikit-learn``, ``scipy``, 
``mygene``, ``fftw``, ``pynndescent``, ``networkx``, ``python-louvain``, ``fitsne``
and a few other standard Python packages. Most of these packages are installed with installation of the 
latest release of ``DigitalCellSorter``:

.. code-block:: bash

    $ pip install DigitalCellSorter

Alternatively, you can install this module directly from GitHub using:

.. code-block:: bash

    $ pip install git+https://github.com/sdomanskyi/DigitalCellSorter

Also one can create a local copy of this project for development purposes by running:

.. code-block:: bash

    $ git clone https://github.com/sdomanskyi/DigitalCellSorter

To install ``fftw`` from the ``conda-forge`` channel add ``conda-forge`` to your channels.
Once the conda-forge channel has been enabled, ``fftw`` can be installed as follows:

.. code-block:: bash 

    $ conda config --add channels conda-forge
    $ conda install fftw

**Loading the package**
-----------------------

In your script import the package:

.. code-block:: python

	import DigitalCellSorter

Create an instance of ``class DigitalCellSorter``. Here, for simplicity, we use Default parameter values:

.. code-block:: python

	DCS = DigitalCellSorter.DigitalCellSorter()
