import os
import urllib.request
import pandas as pd
import DigitalCellSorter

if __name__ == '__main__':
    ans = input('This will execute the demo of DCS. Continue?(y/n): ')
    
    if ans in ['y', 'yes', 'Y']:
        metadata = '''    Data used in this demo:
        Downloaded from https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/manual_5k_pbmc_NGSC3_ch1
        Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Manual (channel 1)
        Single Cell Gene Expression Dataset by Cell Ranger 3.1.0
        PBMCs isolated and cryopreserved by AllCells. PBMCs are primary cells with relatively small amounts of RNA (~1pg RNA/cell).
        Libraries were prepared following the Chromium Single Cell 3ʹ Reagent Kits v3.1 User Guide (CG000204 RevD). 
        8 libraries were prepared on 1 chip. This is channel 1 of the chip.
        5000 cells targeted, 4339 cells detected
        Sequenced on Illumina NovaSeq with approximately 52296 reads per cell
        28bp read1 (16bp Chromium barcode and 12bp UMI), 91bp read2 (transcript), and 8bp I7 sample barcode
        Published on February 28, 2020\n'''

        print(metadata)

        # Get the data. The file will be downloaded from github if not found locally
        try:
            filePath = os.path.join('data', 'manual_5k_pbmc_NGSC3_ch1.gz')

            if not os.path.exists('data'):
                os.makedirs('data')

            if not os.path.isfile(filePath):
                print('Downloading data file')
                urllib.request.urlretrieve('https://github.com/sdomanskyi/DigitalCellSorter/raw/master/data/manual_5k_pbmc_NGSC3_ch1.gz', filePath)
        except Exception as exception:
            print('Could not download the file\n', exception)
            exit()
        else:
            df_data = pd.read_pickle(filePath)

        # Set up a copy of Digital Cell Sorter
        DCS = DigitalCellSorter.DigitalCellSorter(dataName='PBMC5k', geneListFileName='CIBERSORT_LM22_14', verbose=2)
        DCS.saveDir = 'output demo 5k'

        # Process the data, annotate cell types and make plots
        DCS.prepare(df_data)
        DCS.process()
        DCS.annotate()
        DCS.visualize()

        # Do additional analysis and plots
        DCS.makeIndividualGeneTtestPlot('CD4')
        DCS.makeIndividualGeneExpressionPlot(['CD19', DCS.getHugoName('CD56'), DCS.getHugoName('CD16'), 'CD4', 'CD8A', 'CD8B', 'CD33'])
        DCS.makeSankeyDiagram(DCS.getCountsDataframe(DCS.loadAnnotatedLabels(), DCS.loadAnnotatedLabels(infoType='batch')))
        DCS.makeMarkerExpressionPlot(fontscale=0.8)
        DCS.makeHopfieldLandscapePlot(plot3D=True)
        DCS.makeAnomalyScoresPlot()
