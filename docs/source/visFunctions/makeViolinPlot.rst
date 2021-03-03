**Customized violin plot**
==========================

This visualization function can be launched from `class DigitalCellSorter`.
By default the function also exports a summary of the data.

From submodule `VisualizationFunctions`: 


.. autoclass:: DigitalCellSorter.VisualizationFunctions.VisualizationFunctions
    :noindex:

    .. automethod:: makeViolinPlot


**Example with synthetic data:**

.. code-block:: python

    data = pd.DataFrame({'Celltype': np.where((np.random.rand(1000)>0.7), 'Epithelial', 'Endothelial'), 
                         'Condition': np.random.rand(1000)>0.3, 
                         'Gene 1': np.log(np.random.rand(1000)*3 + 1), 
                         'Gene 2': np.log(np.random.rand(1000)*1.5 + 1)})

    DCS = DigitalCellSorter.DigitalCellSorter()
    DCS.makeViolinPlot(data, ['Gene 1', 'Gene 2'], 
                       dimPanels='Celltype', dimCategories='Condition', 
                       title='{name}: {gene}', ylabel='Condition',
                       pointsSize=5, fontsize=12)


**Example output:**

.. thumbnail:: https://github.com/sdomanskyi/DigitalCellSorter/blob/master/docs/examples/violin_plot.png?raw=true
    :title: Example
    :alt: Cannot load this figure
    :align: center
    :width: 250px
    :height: 250px
    :download: false



