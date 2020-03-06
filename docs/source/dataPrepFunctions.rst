**Data preparation functions**
==============================

Set of generic tools for retrieving, loading, and preparation of `Human Cell Atlas (HCA) datasets 
<https://data.humancellatlas.org/explore/projects/>`_ is contained in this module.

Example:

.. code-block:: python
    :emphasize-lines: 1-

    import os
    import DigitalCellSorter.ReadPrepareDataHCA as prep

    # Example URL of a relatvely small dataset of scRNA-seq of human pancreas
    url = "https://data.humancellatlas.org/project-assets/project-matrices/cddab57b-6868-4be4-806f-395ed9dd635a.homo_sapiens.mtx.zip"

    # Path of directories where the data will be placed
    extractPath = os.path.join(os.path.join(os.path.dirname(__file__), ''), 'data', os.path.splitext(os.path.basename(url))[0])

    # Download data and unpack it to a specified directory
    prep.getHCAdataByURL(url, extractPath)

    # Record *.h5 files of individual donor IDs
    IDs = prep.recordFilesOfIndividualDonors(extractPath, organName='islet of Langerhans')

    # Load ready-to-use dataset of the first donor ID
    df = prep.getDataframeByDonorID(extractPath, IDs[0])

    # Print the shape of just loaded dataset
    print(df.shape)


**Submodule ReadPrepareDataHCA**

.. automodule:: DigitalCellSorter.ReadPrepareDataHCA
    :members:
    :undoc-members:
    :show-inheritance:






