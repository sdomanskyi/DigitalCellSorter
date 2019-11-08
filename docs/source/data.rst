**Input Data Format**
=====================

Gene Expression Data Format


The input gene expression data is expected in one of the following formats:

1. Spreadsheet of comma-separated values ``csv`` containing condensed matrix in a form ``('cell', 'gene', 'expr')``. 
If there are batches in the data the matrix has to be of the form ``('batch', 'cell', 'gene', 'expr')``. Columns order can be arbitrary.


+------+------+------+
| cell | gene | expr |
+======+======+======+
| C1   | G1   | 3    |
+------+------+------+
| C1   | G2   | 2    |
+------+------+------+
| C1   | G3   | 1    |
+------+------+------+
| C2   | G1   | 1    |
+------+------+------+
| C2   | G4   | 5    |
+------+------+------+
| ...  | ...  | ...  |
+------+------+------+

or:

+--------+------+------+------+
| batch  | cell | gene | expr |
+========+======+======+======+
| batch0 | C1   | G1   | 3    |
+--------+------+------+------+
| batch0 | C1   | G2   | 2    |
+--------+------+------+------+
| batch0 | C1   | G3   | 1    |
+--------+------+------+------+
| batch1 | C2   | G1   | 1    |
+--------+------+------+------+
| batch1 | C2   | G4   | 5    |
+--------+------+------+------+
| ...    | ...  | ...  | ...  |
+--------+------+------+------+




2. Spreadsheet of comma-separated values ``csv`` where rows are genes, columns are cells with gene expression counts.
If there are batches in the data the spreadsheet the first row should be ``'batch'`` and the second ``'cell'``.


+-------+--------+--------+--------+--------+
| cell  | C1     | C2     | C3     | C4     |
+=======+========+========+========+========+
| G1    |        | 3      | 1      | 7      |
+-------+--------+--------+--------+--------+
| G2    | 2      | 2      |        | 2      |
+-------+--------+--------+--------+--------+
| G3    | 3      | 1      |        | 5      |
+-------+--------+--------+--------+--------+
| G4    | 10     |        | 5      | 4      |
+-------+--------+--------+--------+--------+
| ...   | ...    | ...    | ...    | ...    |
+-------+--------+--------+--------+--------+

or:

+-------+--------+--------+--------+--------+
| batch | batch0 | batch0 | batch1 | batch1 |
+=======+========+========+========+========+
| cell  | C1     | C2     | C3     | C4     |
+-------+--------+--------+--------+--------+
| G1    |        | 3      | 1      | 7      |
+-------+--------+--------+--------+--------+
| G2    | 2      | 2      |        | 2      |
+-------+--------+--------+--------+--------+
| G3    | 3      | 1      |        | 5      |
+-------+--------+--------+--------+--------+
| G4    | 10     |        | 5      | 4      |
+-------+--------+--------+--------+--------+
| ...   | ...    | ...    | ...    | ...    |
+-------+--------+--------+--------+--------+

3. ``Pandas DataFrame`` where ``axis 0`` is genes and ``axis 1`` are cells.
If the are batched in the data then the index of ``axis 1`` should have two levels, e.g. ``('batch', 'cell')``, 
with the first level indicating patient, batch or expreriment where that cell was sequenced, and the
second level containing cell barcodes for identification.

.. code:: python

    df = pd.DataFrame(data=[[2,np.nan],[3,8],[3,5],[np.nan,1]], 
                      index=['G1','G2','G3','G4'], 
                      columns=pd.MultiIndex.from_arrays([['batch0','batch1'],['C1','C2']], names=['batch', 'cell']))    


4. ``Pandas Series`` where index should have two levels, e.g. ``('cell', 'gene')``. If there are batched in the data
the first level should be indicating patient, batch or expreriment where that cell was sequenced, the second level cell barcodes for 
identification and the third level gene names.

.. code:: python

    se = pd.Series(data=[1,8,3,5,5], 
                   index=pd.MultiIndex.from_arrays([['batch0','batch0','batch1','batch1','batch1'],
                                                    ['C1','C1','C1','C2','C2'],
                                                    ['G1','G2','G3','G1','G4']], names=['batch', 'cell', 'gene']))


Any of the data types outlined above need to be prepared/validated with a function ``prepare()``. 

