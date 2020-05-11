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

For detailed instructions and other ways to install ``DigitalCellSorter`` as wells as
list of optional packages and instructions on how to install them see
**Prerequisites** section at https://github.com/sdomanskyi/DigitalCellSorter


**Loading the package**
-----------------------

In your script import the package:

.. code-block:: python

	import DigitalCellSorter

Create an instance of ``class DigitalCellSorter``. Here, for simplicity, we use Default parameter values:

.. code-block:: python

	DCS = DigitalCellSorter.DigitalCellSorter()
