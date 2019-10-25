import pandas as pd
import numpy as np
import DigitalCellSorter.DigitalCellSorter as DigitalCellSorter
    
if __name__ == '__main__':

    DCS = DigitalCellSorter.DigitalCellSorter()

    # PBMC DCS TD6 vs DropClust 
    if True:
        with open('ColormapForCellTypes_TD_coarse_6.txt', 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}
        colormap['Missing'] = (0.0, 0.0, 0.0, 1.0)
        colormap = DCS.convertColormap(colormap)

        se_Drop = pd.read_excel('predictions.xlsx', sheet_name='predictions', header=0, index_col=0)['annotation by dropClust']
        se_DCS = pd.read_hdf('PBMC_processed_TD_ext_6.h5', key='df_markers_expr', mode='r').T.reset_index().set_index('cell')['label'].str.split(' #', expand=True)[0]
    
        df_new_marker_list = DCS.getCountsDataframe(se_DCS, se_Drop)

        typesMapping1 = {item:item.split(' #')[0] for item in np.unique(se_DCS)}
        typesMapping1.update({'Missing':'Missing'})

        typesMapping2 = {'Naive T cells':'CD.4.T.cells', 
                        'CD4+ Memory T cells':'CD.4.T.cells', 
                        'Natural Killer T (NKT) cells':'NK.cells', 
                        'B cells':'B.cells', 
                        'CD8+ T cells':'CD.8.T.cells', 
                        'Natural Killers':'NK.cells', 
                        'CD8+ T cells':'CD.8.T.cells',
                        'CD16+ Monocytes':'monocytic.lineage', 
                        'CD14+ Monocytes':'monocytic.lineage', 
                        'Regulatory T (Treg) cells':'CD.4.T.cells', 
                        'Monocyte Derived Dendritic Cells':'Unknown', 
                        'Circulating Megakaryocyte Progenitors':'Unknown', 
                        'Natural Killer Progenitors (NKP)':'NK.cells', 
                        'Plasmacytoid Dendritic Cells':'monocytic.lineage',
                        'Missing':'Missing'}

        colormapForIndex = {k:colormap[v] for k, v in typesMapping1.items()}
        colormapForColumns = {k:colormap[v] for k, v in typesMapping2.items()}

        DCS.makeSankeyDiagram(df_new_marker_list, 'DCS_vs_DropClust ext 6', colormapForIndex, colormapForColumns, title='DCS TD6 vs DropClust', quality=6, interactive=False)

    # PBMC DCS vs DropClust
    if False:
        with open('ColormapForCellTypes.txt', 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}
        colormap['Missing'] = (0.0, 0.0, 0.0, 1.0)
        colormap = DCS.convertColormap(colormap)

        ##se_DCS = pd.read_hdf('PBMC_processed.h5', key='labels', mode='r').reset_index().set_index('cell')['label']
        se_DCS = pd.read_hdf('PBMC_processed.h5', key='df_markers_expr', mode='r').T.reset_index().set_index('cell')['label'].str.split(' #', expand=True)[0]
        #se_DCS = pd.read_hdf('PBMC_processed.h5', key='labels', mode='r').reset_index().set_index('cell')['label'].str.split(' #', expand=True)[0]
        se_Drop = pd.read_excel('predictions by dropCLust.xlsx', header=None, index_col=0)[1]
    
        df_new_marker_list = DCS.getCountsDataframe(se_DCS, se_Drop)

        typesMapping1 = {item:item.split(' #')[0] for item in np.unique(se_DCS)}
        typesMapping2 = {'B cell': 'B cell', 
                 'Cytotoxic T cell': 'T cell', 
                 'Dendritic cell': 'Dendritic cell', 
                 'Megakaryocyte': 'Unknown', 
                 'Memory T cell': 'T cell', 
                 'Monocyte': 'Monocytes/Macrophages', 
                 'NK T cell': 'NK cell', 
                 'NK cell': 'NK cell', 
                 'Naive T cell':  'T cell', 
                 'Regulatory T (Treg) cell': 'T cell'}

        typesMapping1.update({'Missing':'Missing'})
        typesMapping2.update({'Missing':'Missing'})

        colormapForIndex = {k:colormap[v] for k, v in typesMapping1.items()}
        colormapForColumns = {k:colormap[v] for k, v in typesMapping2.items()}

        DCS.makeSankeyDiagram(df_new_marker_list, 'DCS_vs_DropClust detailed 2', colormapForIndex, colormapForColumns, title='DCS vs DropClust', quality=6, interactive=True)

    # BM1 clustering KMeans vs Agglomerative
    if False:
        with open('ColormapForCellTypes.txt', 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}
        colormap['Missing'] = (0.0, 0.0, 0.0, 1.0)
        colormap = DCS.convertColormap(colormap)

        se1 = pd.read_hdf('BM1_processed_KMeans.h5', key='df_markers_expr', mode='r').columns.to_series().reset_index().set_index('cell')['label'].str.split(' #', expand=True)[0]
        #se_DCS_1 = pd.read_hdf('BM1_processed_KMeans.h5', key='df_markers_expr', mode='r').columns.to_series().reset_index().set_index('cell')['label']
        se2 = pd.read_hdf('BM1_processed_Aggl.h5', key='df_markers_expr', mode='r').columns.to_series().reset_index().set_index('cell')['label'].str.split(' #', expand=True)[0]
        #se_DCS_2 = pd.read_hdf('BM1_processed_Aggl.h5', key='df_markers_expr', mode='r').columns.to_series().reset_index().set_index('cell')['label']
    
        df_new_marker_list = DCS.getCountsDataframe(se1, se2)

        typesMapping1 = {item:item.split(' #')[0] for item in np.unique(se1)}
        typesMapping2 = {item:item.split(' #')[0] for item in np.unique(se2)}
        typesMapping1.update({'Missing':'Missing'})
        typesMapping2.update({'Missing':'Missing'})

        colormapForIndex = {k:colormap[v] for k, v in typesMapping1.items()}
        colormapForColumns = {k:colormap[v] for k, v in typesMapping2.items()}

        DCS.makeSankeyDiagram(df_new_marker_list, 'BM1Kmeans_vs_BM1Aggl detailed', colormapForIndex, colormapForColumns, title='BM1: Kmeans clustering vs Agglomerative clustering', quality=6)

    # BM1 gene list CIBERSORT vs own
    if False:
        name1 = 'BM1_processed_Aggl'
        name2 = 'BM1_processed_top50'
        fileName = 'BM1 CIBERSORT_vs_our merged'
        title = 'BM1: CIBERSORT gene list vs our'
        quality = 6
        detailed = False

        with open('ColormapForCellTypes.txt', 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}
        colormap['Missing'] = (0.0, 0.0, 0.0, 1.0)
        colormap = DCS.convertColormap(colormap)

        se1 = pd.read_hdf(name1 + '.h5', key='df_markers_expr', mode='r').columns.to_series().reset_index().set_index('cell')['label']
        se2 = pd.read_hdf(name2 + '.h5', key='df_markers_expr', mode='r').columns.to_series().reset_index().set_index('cell')['label']

        if not detailed:
            se1 = se1.str.split(' #', expand=True)[0]
            se2 = se2.str.split(' #', expand=True)[0]
    
        df_new_marker_list = DCS.getCountsDataframe(se1, se2)

        typesMapping1 = {item:item.split(' #')[0] for item in np.unique(se1)}
        typesMapping1.update({'Missing':'Missing'})
        colormapForIndex = {k:colormap[v] for k, v in typesMapping1.items()}

        typesMapping2 = {item:item.split(' #')[0] for item in np.unique(se2)}
        connection = {'B.cells':'B cell',
                        'CD4.T.cells':'T cell',
                        'CD8.T.cells':'T cell',
                        'NK.cells':'NK cell',
                        'neutrophils':'Granulocyte',
                        'monocytic.lineage':'Monocytes/Macrophages',
                        'fibroblasts':'Dendritic cell',
                        'endothelial.cells':'Dendritic cell',
                        'Unknown':'Unknown',
                        'Missing':'Missing'}
        typesMapping2.update({'Missing':'Missing'})
        colormapForColumns = {k:colormap[connection[v]] for k, v in typesMapping2.items()}

        DCS.makeSankeyDiagram(df_new_marker_list, fileName + (' detailed' if detailed else ' merged'), colormapForIndex, colormapForColumns, title=title, quality=quality)

    # BM123 COMBAT
    if False:
        name1 = 'BM123_no_corr_processed.h5'
        name2 = 'BM123_with_corr_processed.h5'
        fileName = 'no corr vs corr'
        title = 'BM: no corr vs corr'
        quality = 6
        detailed = True

        with open('ColormapForCellTypes.txt', 'r') as temp_file:
            colormap = {item.strip().split('\t')[0]:eval(item.strip().split('\t')[1]) for item in temp_file.readlines()}
        colormap['Missing'] = (0.0, 0.0, 0.0, 1.0)
        colormap = DCS.convertColormap(colormap)

        se1 = pd.read_hdf(name1, key='df_markers_expr', mode='r').columns.to_series().reset_index().set_index('cell')['label']
        se2 = pd.read_hdf(name2, key='df_markers_expr', mode='r').columns.to_series().reset_index().set_index('cell')['label']

        if not detailed:
            se1 = se1.str.split(' #', expand=True)[0]
            se2 = se2.str.split(' #', expand=True)[0]
    
        df_new_marker_list = DCS.getCountsDataframe(se1, se2)

        connection = {'B.cells':'B cell',
                        'CD4.T.cells':'T cell',
                        'CD8.T.cells':'T cell',
                        'NK.cells':'NK cell',
                        'neutrophils':'Granulocyte',
                        'monocytic.lineage':'Monocytes/Macrophages',
                        'fibroblasts':'Dendritic cell',
                        'endothelial.cells':'Dendritic cell',
                        'Unknown':'Unknown',
                        'Missing':'Missing'}

        typesMapping1 = {item:item.split(' #')[0] for item in np.unique(se1)}
        typesMapping1.update({'Missing':'Missing'})
        colormapForIndex = {k:colormap[v] for k, v in typesMapping1.items()}

        typesMapping2 = {item:item.split(' #')[0] for item in np.unique(se2)}
        typesMapping2.update({'Missing':'Missing'})
        colormapForColumns = {k:colormap[v] for k, v in typesMapping2.items()}

        DCS.makeSankeyDiagram(df_new_marker_list, fileName + (' detailed' if detailed else ' merged'), colormapForIndex, colormapForColumns, title=title, quality=quality)



        
    ## PBMC ################################################################################################################################################
    #if False:
    #    dataName = 'filtered_matrices_mex'
    #    df_expr = pd.read_csv(os.path.join('data', dataName, 'matrix.mtx'), header=None, skiprows=3)[0].str.strip().str.split(' ', expand=True).applymap(int).set_index([0,1])
    #    df_expr.index.names = ['gene', 'cell']
    #    df_expr.columns = ['PBMC']
    #    df_expr = df_expr.unstack(fill_value=0)
    #    df_expr.columns.names = ['batch', 'cell']
    #    df_expr.index = pd.read_csv(os.path.join('data', dataName, 'genes.tsv'), delimiter='\t', header=None)[1].loc[df_expr.index-1]
    #    cells = pd.read_csv(os.path.join('data', dataName, 'barcodes.tsv'), delimiter='\t', header=None)[0].loc[df_expr.columns.get_level_values('cell')-1]
    #    df_expr.columns = pd.MultiIndex.from_arrays([df_expr.columns.get_level_values('batch'), cells], names=['batch', 'cell'])
    #    df_expr.to_hdf(os.path.join('data', 'PBMC.h5'), key='PBMC', mode='a', complevel=4, complib='zlib')
    #else:
    #    dataName = 'PBMC'
    #    df_expr = pd.read_hdf(os.path.join('data', 'PBMC.h5'), key='PBMC', mode='r')

    ## BM ##################################################################################################################################################
    #df_expr = pd.DataFrame()
    #for dataName in ['BM1', 'BM2', 'BM3']:
    #    #HCA.PrepareDataOnePatient(os.path.join('data', 'ica_bone_marrow_h5.h5'), dataName, os.path.join('data', ''), useAllData=False if os.name == 'nt' else True, cellsLimitToUse=1000)
    #    df_expr = pd.concat([df_expr, pd.read_hdf(os.path.join('data', 'HCA_%s.h5'%(dataName)), key=dataName, mode='r')], sort=False, axis=1)
    #    DigitalCellSorter.timeMark()