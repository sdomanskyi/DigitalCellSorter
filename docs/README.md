# Documentation of Digital Cell Sorter Toolkit

[![DOI](https://badge.fury.io/gh/sdomanskyi%2FDigitalCellSorter.svg)](https://badge.fury.io/gh/sdomanskyi%2FDigitalCellSorter)
[![DOI](https://badge.fury.io/py/DigitalCellSorter.svg)](https://pypi.org/project/DigitalCellSorter)
[![DOI](https://readthedocs.org/projects/digital-cell-sorter/badge/?version=latest)](https://digital-cell-sorter.readthedocs.io/en/latest/?badge=latest)

The ready-to-read documentation is available at https://digital-cell-sorter.readthedocs.io/.
The documentation of our software is built with [**Sphinx**](https://www.sphinx-doc.org/ "Sphinx") at 
[**ReadTheDocs.org**](https://readthedocs.org/).

## Documentation Updates

Any changes to source code and docstrings
are automatically reflected at [**ReadTheDocs.org**](https://readthedocs.org/) 
(a new version of documentation is built silently). 

## Documentation Development

For development and testing of the documentation locally (on the development machine) 
install [**Sphinx**](https://www.sphinx-doc.org/ "Sphinx") by:

	$ pip install -U sphinx

To compile **html** version of the documentation go to **docs/** directory and run:

	$ sphinx-build -E -a -b html ./source ./build

We are utilizing a 3rd party **Sphinx** extension [**sphinxcontrib-images**](https://github.com/sphinx-contrib/images 
"GitHub repository") extension, allowing to display documentation images in a organized way.


