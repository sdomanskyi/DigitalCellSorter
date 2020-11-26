**Data preparation**
==============================


**Output from kallisto-bustools (kp-python)**
---------------------------------------------

In this example we use raw sequencing data stored in ``.fastq`` format, from 1000 PBMC, 
the data can be accessed at
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3

.. note:: This is by no means a tutorial for processing scRNA-seq data. We only demonstrate the workflow of connecting upstream analysis software and DCS.

Download the data and unpack the ``.tar`` file (~5.17 GB):

::

   wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
   tar -xvf pbmc_1k_v3_fastqs.tar

To process sequencing data one could use ``kallisto bus`` tool to generate BUS file
following by ``bustools count`` to generate count matrices from a BUS file. 
However, we prefer to use ``kb-python``, a package that wraps the ``kallisto``
and ``bustools`` single-cell RNA-seq workflow (Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature biotechnology, 34(5), 525)
``kb-python`` can be installed with ``pip``.

::

    kb count -i kallisto_index/homo_sapiens/transcriptome.idx \
             -g kallisto_index/homo_sapiens/transcripts_to_genes.txt \
             -x 10xv3 \
             --filter \
             -t 4 \
             pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz \
             pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz \
             pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz \
             pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz

.. container:: toggle

    .. container:: header

       **Output from kb count command above**
    [2020-11-20 14:32:51,136]    INFO Using index kallisto_index/homo_sapiens/transcriptome.idx to generate BUS file to . from
    
    [2020-11-20 14:32:51,136]    INFO         pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz
    
    [2020-11-20 14:32:51,136]    INFO         pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz
    
    [2020-11-20 14:32:51,136]    INFO         pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz
    
    [2020-11-20 14:32:51,136]    INFO         pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz
    
    [2020-11-20 14:36:33,477]    INFO Sorting BUS file ./output.bus to ./tmp/output.s.bus
    
    [2020-11-20 14:36:57,118]    INFO Whitelist not provided
    
    [2020-11-20 14:36:57,118]    INFO Copying pre-packaged 10XV3 whitelist to .
    
    [2020-11-20 14:36:57,675]    INFO Inspecting BUS file ./tmp/output.s.bus
    
    [2020-11-20 14:37:07,641]    INFO Correcting BUS records in ./tmp/output.s.bus to ./tmp/output.s.c.bus with whitelist ./10xv3_whitelist.txt
    
    [2020-11-20 14:37:29,264]    INFO Sorting BUS file ./tmp/output.s.c.bus to ./output.unfiltered.bus
    
    [2020-11-20 14:37:47,478]    INFO Generating count matrix ./counts_unfiltered/cells_x_genes from BUS file ./output.unfiltered.bus
    
    [2020-11-20 14:37:59,662]    INFO Filtering with bustools
    
    [2020-11-20 14:37:59,662]    INFO Generating whitelist ./filter_barcodes.txt from BUS file ./output.unfiltered.bus
    
    [2020-11-20 14:37:59,790]    INFO Correcting BUS records in ./output.unfiltered.bus to ./tmp/output.unfiltered.c.bus with whitelist ./filter_barcodes.txt
    
    [2020-11-20 14:38:14,344]    INFO Sorting BUS file ./tmp/output.unfiltered.c.bus to ./output.filtered.bus
    
    [2020-11-20 14:38:30,918]    INFO Generating count matrix ./counts_filtered/cells_x_genes from BUS file ./output.filtered.bus

|

The output directory that we are interested in is ``counts_filtered/``. Rename it:

::

   mv counts_filtered/ kb_1k_PBMC_output/

This will generate counts data in the directory ``kb_1k_PBMC_output/``.


**Output from CellRanger**
--------------------------

Here we use CellRanger-processed data stored in ``.mtx`` format, from 1000 PBMC, 
the data can be accessed at
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3

Download the data and unpack the ``.tar.gz`` file (~9 MB):

::

   wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz
   tar -xzf pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz && mv filtered_feature_bc_matrix/ cellranger_1k_PBMC_output/

These two commands will prepare the processed counts data in the directory ``cellranger_1k_PBMC_output/``.



**Import from kallisto-bustools (kp-python)**
---------------------------------------------

.. code-block:: python
    :emphasize-lines: 5

    import DigitalCellSorter
    from DigitalCellSorter.core import readMTXdata

    # Read the MTX data
    df = readMTXdata(dataDir='kb_1k_PBMC_output/', origin='kb-python')

    # (Optional) Convert gene names to HUGO
    DCS = DigitalCellSorter.DigitalCellSorter()
    DCS.prepare(df)
    DCS.convert('ensembl', 'hugo')

    # Check the DCS data
    print(DCS.df_expr)


**Import from CellRanger**
--------------------------

.. code-block:: python
    :emphasize-lines: 5

    import DigitalCellSorter
    from DigitalCellSorter.core import readMTXdata

    # Read the MTX data
    df = readMTXdata(dataDir='cellranger_1k_PBMC_output/', origin='cellranger')        

    # (Optional) Convert gene names to HUGO
    DCS = DigitalCellSorter.DigitalCellSorter()
    DCS.prepare(df)
    DCS.convert('ensembl', 'hugo')

    # Check the DCS data
    print(DCS.df_expr)



**Function readMTXdata**
--------------------------

Function to read data in MTX format (see usage examples above).

.. autofunction:: DigitalCellSorter.core.readMTXdata



**Human Cell Atlas tools**
--------------------------

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






