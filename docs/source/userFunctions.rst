**User Functions**
==================

User functions from **DigitalCellSorter.core.DigitalCellSorter** class.

.. Note:: All of the tools listed below in this section are intended to use from an 
    instance of a ``DigitalCellSorter`` class. For example:

    .. code-block:: python
        :emphasize-lines: 1-

        DCS = DigitalCellSorter.DigitalCellSorter()

        DCS.dataName = 'my_data_name'
        DCS.saveDir = os.path.join(os.path.dirname(__file__), 'output', DCS.dataName, '')

        data = DCS.prepare(raw_data)

        DCS.process(DCS.prepare(data))

        DCS.getIndividualGeneExpressionPlot('CCL5')

        DCS.getIndividualGeneTtestPlot('CCL5', analyzeBy='celltype')

        cells = DCS.getCells(celltype='T cell')
        DCS.makeAnomalyScoresPlot(cells)

        # ...

    Direct use of function from where they are stored may result in 
    undefined behavior.

.. container:: toggle

    .. container:: header

       **Description of the package functionality**

    .. automodule:: DigitalCellSorter.core
        :exclude-members: DigitalCellSorter
        :noindex:

|
|

**Primary tools**
-----------------

Primary tools are used for pre-processing of the input data, quality control, batch 
correction, dimensionality reduction, clustering and cell type annotation.

.. Note:: We reccomend to use only functions ``prepare()``, ``process()``, and 
     ``visualize()`` of the Primary tools. All processing workflow is contained 
     within ``process()``.
     If you wish to modify the workflow use the other components of the 
     Primary tools, such as ``cluster()``, ``project()`` etc. 

.. currentmodule:: DigitalCellSorter.core.DigitalCellSorter

**References to DigitalCellSorter class:**

.. autosummary:: 
     prepare
     convert
     clean
     project
     cluster
     vote
     process
     visualize


.. container:: toggle

    .. container:: header

       **Function** ``prepare()``: prepare input data for function ``process()``

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.prepare
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``convert()``: convert gene index of a DataFrame prepared by function ``prepare()``
       from one naming convention to another

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.convert
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``clean()``: validate index, replace missing with zeros, 
       remove all-zero rows and columns of a DataFrame

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.clean
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``normalize()``: rescale all cells, log-transform data,
       remove constant genes, and sort index of a DataFrame

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.normalize
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``project()``: project data to lower dimensions

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.project
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``cluster()``: cluster PCA-reduced data into a desired number of clusters

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.cluster
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``vote()``: produce cluster voting results

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.vote
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``process()``: main function

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.process
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``visualize()``: make all default plots of to visualize results 
       of function ``process()``

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.visualize
     	:noindex:

|


**Extraction tools**
--------------------

.. warning:: Use these functions only after ``process()``


**References to DigitalCellSorter class:**

.. currentmodule:: DigitalCellSorter.core.DigitalCellSorter
.. autosummary:: 
     getExprOfGene
     getExprOfCells
     getCells
     getAnomalyScores
     getNewMarkerGenes
     getIndexOfGoodQualityCells
     getCountsDataframe


.. container:: toggle

    .. container:: header

       **Function** ``getExprOfGene()``: Get expression of a gene

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getExprOfGene
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``getExprOfCells()``: Get expression of a set of cells

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getExprOfCells
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``getCells()``: get cells index by celltype, clusterIndex or clusterName

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getCells
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``getAnomalyScores()``: get anomaly score of cells based on some reference set

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getAnomalyScores
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``getNewMarkerGenes()``: extract new markers from the annotated clusters and produce plot of the new markers 

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getNewMarkerGenes
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``getIndexOfGoodQualityCells()``: Get index of sells that satisfy the QC criteria

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getIndexOfGoodQualityCells
     	:noindex:

|

.. container:: toggle

    .. container:: header

       **Function** ``getCountsDataframe()``: Get a pandas.DataFrame with cross-counts (overlaps) between two pandas.Series

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getCountsDataframe
     	:noindex:

|



**Visualization tools**
-----------------------

.. warning:: Use these functions only after ``process()``


**References to DigitalCellSorter class:**

.. currentmodule:: DigitalCellSorter.core.DigitalCellSorter
.. autosummary:: 
     getProjectionPlotAnnotated
     getProjectionPlotByBatches
     getProjectionPlotByClusters
     getProjectionPlotsQC
     getMarkerSubplots
     getAnomalyScoresPlot
     getIndividualGeneTtestPlot
     getIndividualGeneExpressionPlot

**References to VisualizationFunctions class:**

.. currentmodule::  DigitalCellSorter.VisualizationFunctions.VisualizationFunctions
.. autosummary:: 
     makeQualityControlHistogramPlot
     makeHistogramNullDistributionPlot
     makeVotingResultsMatrixPlot
     makeMarkerExpressionPlot
     makeStackedBarplot
     makeSankeyDiagram


.. container:: toggle

    .. container:: header

       **Function** ``getProjectionPlotAnnotated()``: Produce t-SNE plot colored by cell types

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getProjectionPlotAnnotated
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_clusters_annotated.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 250px
    :height: 250px
    :download: false

|

.. container:: toggle

    .. container:: header

       **Function** ``getProjectionPlotByBatches()``: Produce t-SNE plot colored by batches

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getProjectionPlotByBatches
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_patients.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 250px
    :height: 250px
    :download: false

|

.. container:: toggle

    .. container:: header

       **Function** ``getProjectionPlotByClusters()``: Produce t-SNE plot colored by clusters

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getProjectionPlotByClusters
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_clusters.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 250px
    :height: 250px
    :download: false

|

.. container:: toggle

    .. container:: header

       **Function** ``makeQualityControlHistogramPlot()``: Produce Quality Control histogram plots

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.makeQualityControlHistogramPlot
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/QC_plots/BM1_number_of_genes_histogram.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 200px
    :height: 200px
    :download: false
    :group: QC plots

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/QC_plots/BM1_count_depth_histogram.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 200px
    :height: 200px
    :download: false
    :group: QC plots

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/QC_plots/BM1_fraction_of_mitochondrialGenes_histogram.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 200px
    :height: 200px
    :download: false
    :group: QC plots

.. container:: clearfix

   .. stuff

|

.. container:: toggle

    .. container:: header

       **Function** ``getProjectionPlotsQC()``: Produce Quality Control t-SNE plots

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getProjectionPlotsQC
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_number_of_genes.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 200px
    :height: 200px
    :download: false
    :group: QC plots

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_count_depth.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 200px
    :height: 200px
    :download: false
    :group: QC plots

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_fraction_of_mitochondrialGenes.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 200px
    :height: 200px
    :download: false
    :group: QC plots

.. container:: clearfix

   .. stuff

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_is_quality_cell.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 300px
    :height: 300px
    :download: false
    :group: QC plots

|

.. container:: toggle

    .. container:: header

       **Function** ``getMarkerSubplots()``: Produce subplots on each marker and its expression on all clusters

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getMarkerSubplots
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/marker_subplots/BM1_CD4_(CD4_CD4mut).png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 300px
    :height: 300px
    :download: false
    :group: marker subplots

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/marker_subplots/BM1_CD19_(B4_CVID3_CD19).png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 300px
    :height: 300px
    :download: false
    :group: marker subplots

.. container:: clearfix

   .. stuff

|

.. container:: toggle

    .. container:: header

       **Function** ``getAnomalyScoresPlot()``: Make anomaly scores plot

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getAnomalyScoresPlot
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_anomaly_score%20All.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 300px
    :height: 300px
    :download: false
    :group: anomaly

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_anomaly_score%20Cluster2.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 300px
    :height: 300px
    :download: false
    :group: anomaly

.. container:: clearfix

   .. stuff

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_anomaly_score%20B%20cell.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 300px
    :height: 300px
    :download: false
    :group: anomaly

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_clusters_by_anomaly_score%20T%20cell.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: left
    :width: 300px
    :height: 300px
    :download: false
    :group: anomaly

.. container:: clearfix

   .. stuff

|

.. container:: toggle

    .. container:: header

       **Function** ``getIndividualGeneTtestPlot()``: Produce individual gene t-test plot of the two-tailed p-value

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getIndividualGeneTtestPlot
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_ttest_CD4_(CD4_CD4mut).png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 250px
    :height: 250px
    :download: false

|

.. container:: toggle

    .. container:: header

       **Function** ``getIndividualGeneExpressionPlot()``: Produce individual gene expression plot on a 2D layout

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.getIndividualGeneExpressionPlot
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/marker_subplots/BM1_CD4_(CD4_CD4mut).png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 250px
    :height: 250px
    :download: false

|

.. container:: toggle

    .. container:: header

       **Function** ``makeHistogramNullDistributionPlot()``: Produce histogram plot of the voting null distributions

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.makeHistogramNullDistributionPlot
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_null_distributions.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 400px
    :height: 400px
    :download: false

|

.. container:: toggle

    .. container:: header

       **Function** ``makeVotingResultsMatrixPlot()``: Produce voting results voting matrix plot

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.makeVotingResultsMatrixPlot
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_matrix_voting.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 400px
    :height: 400px
    :download: false

|

.. container:: toggle

    .. container:: header

       **Function** ``makeMarkerExpressionPlot()``: Produce image on marker genes and their expression on all clusters

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.makeMarkerExpressionPlot
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_voting.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 600px
    :height: 250px
    :download: false

|

.. container:: toggle

    .. container:: header

       **Function** ``makeStackedBarplot()``: Produce stacked barplot with cell fractions

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.makeStackedBarplot
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/output/BM1/BM1_subclustering_stacked_barplot_.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 150px
    :height: 250px
    :download: false

|

.. container:: toggle

    .. container:: header

       **Function** ``makeSankeyDiagram()``: Make a Sankey diagram, also known as 'river plot' with two groups of nodes

    .. automethod:: DigitalCellSorter.core.DigitalCellSorter.makeSankeyDiagram
     	:noindex:

Example output:

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/Sankey_example.png?raw=true
    :title: Example
    :alt: Cannot load this photo
    :align: center
    :width: 400px
    :height: 400px
    :download: false

|





