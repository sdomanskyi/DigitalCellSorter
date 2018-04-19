# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""



from Bio import Entrez
import time
import mygene
import textwrap

def GetGeneSummary(geneIDs, #iterable containing Hugo and/or Entrez IDs
                   nicePrint = False, #if true, print a nicely formatted summary for each gene
                   returnDict = True, #if true, return a dictionary with gene IDs as the keys and summaries as the values
                   ):
    singleInputFlag = False
    if type(geneIDs) is str: geneIDs = [geneIDs]; singleInputFlag = True
    mg = mygene.MyGeneInfo()
    ezIDs = [mg.query(geneID,size=1,species='human')['hits'][0]['_id'] for geneID in geneIDs]
    hugoIDs = [mg.query(ezID,size=1,species='human')['hits'][0]['symbol'] for ezID in ezIDs]
    Entrez.email='szedlak1@msu.edu'
    request = Entrez.epost('gene',id=','.join(ezIDs),)
    result = Entrez.read(request)
    data = Entrez.esummary(db="gene", webenv=result["WebEnv"], query_key=result["QueryKey"])
    failCounter = 0
    while True:
        try:
            annotations = Entrez.read(data)
            break
        except RuntimeError:
            failCounter+=1
            if failCounter==3:
                print 'Failed to connect! Something''s wrong.'
                break
            print 'Failed to connect (maybe exceeded rate limit). Sleeping and retrying...'
            time.sleep(3)
    summaries = [annotations['DocumentSummarySet']['DocumentSummary'][i]['Summary'] for i in xrange(len(ezIDs))]
    summaries = [summary if len(summary)>0 else 'NCBI summary empty.' for summary in summaries]
    summaries = [summary + ' URL: https://www.ncbi.nlm.nih.gov/gene/' + str(ezIDs[i]) + ' ' for i,summary in enumerate(summaries)]
    if nicePrint:
        summaries_nice = [textwrap.fill(summary, 100) for summary in summaries]
        for geneID,hugoID,summary_nice,ezID in zip(geneIDs,hugoIDs,summaries_nice,ezIDs):
            try: int(geneID); flag = True
            except ValueError: flag = False
            if geneID==hugoID or flag: print '<<< %s >>> (EzID: %s)\n%s\n'%(hugoID,ezID,summary_nice)
            else: print '<<< %s (commonly %s) >>> (EzID: %s)\n%s\n'%(geneID,hugoID,ezID,summary_nice)
    if singleInputFlag==False:
        if returnDict:
            summaryDict = {}
            summaryDict.update(zip([int(ezID) for ezID in ezIDs],summaries))
            summaryDict.update(zip(hugoIDs,summaries))
            summaryDict.update(zip(geneIDs,summaries))
            return summaryDict
        else: return summaries
    else: return summaries[0]



if __name__=='__main__':
    #Define list of gene IDs (can be mixed list of HugoIDs and EzIDs). Can also
    #give single gene as input (e.g. geneIDs='ABL1'), but you should *always*
    #try to do everything at once because the server has rate limits. This is 
    #<<BAD>> practice:
    #
    #  summaryDict = {}
    #  for geneID in ['ABL1','FGFR4',1002]:
    #      summaryDict[geneID] = GetGeneSummary(geneID)
    #
    #The following is good practice:
    geneIDs = ['ABL1',1002,'MNK2']
    #Get a dictionary with the NCBI summary for each gene in the list
    summaryDict = GetGeneSummary(geneIDs,nicePrint=True)
    #The dictionary takes either the HugoID or EzID as the input key 
    #(regardless of which convention was used in the initial query) and returns
    #the NCBI summary 
    print ''
    print '====================='
    print ''
    print 'EzID 25 (ABL1): '+summaryDict[25]






















