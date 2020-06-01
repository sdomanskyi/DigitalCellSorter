import os
import pandas as pd
import numpy as np
from .GenericFunctions import read, write
import mygene

class GeneNameConverter:
    
    def __init__(self, dictDir = None, jumpStart = True, species = 'Human'):

        self.species = species

        if self.species == 'Human':
            self.ensemblStr = 'ENSG' 
        elif self.species == 'Mouse':
            self.ensemblStr = 'ENSM' 
        else:
            raise NotImplementedError

        if dictDir is None: 
            self.dictDir = os.path.join('ensembl_hugo_entrez_alias_dict' % (self.species))
        else: 
            self.dictDir = os.path.join(dictDir, 'ensembl_hugo_entrez_alias_dict_%s' % (self.species))

        #Try getting dictionary from disk, and if it can't be found, initialize an empty dict
        try:
            self.conversionDict = read(self.dictDir)
        except Exception as exception:
            print(exception)

            self.conversionDict = {'hugo': {'entrez': {}, 'ensembl': {}, 'alias': {}},
                                   'entrez': {'hugo': {}, 'ensembl': {}, 'retired': {}},
                                   'ensembl': {'entrez': {}, 'hugo': {}},
                                   'alias': {'hugo': {}, 'entrez': {}, 'ensembl': {}},
                                   'retired': {'entrez': {}}
                                   }
            
            for sourceType in ['hugo', 'alias', 'entrez', 'ensembl', 'retired']:
                self.conversionDict[sourceType]['known'] = set()
                self.conversionDict[sourceType]['unknown'] = set()

            if jumpStart:
                try: 
                    self.JumpStart()
                except IOError: 
                    pass
            else:
                write(self.conversionDict, self.dictDir)

        return
    
    def JumpStart(self,):
        
        #Get a sandwich because this will take a while
        ensemblList = np.loadtxt(os.path.join(os.path.dirname(self.dictDir), 'geneLists', '%s_ensembl.gz') % (self.species), dtype=str, delimiter='\t').flatten()[1:].tolist()
        symbolList = np.loadtxt(os.path.join(os.path.dirname(self.dictDir), 'geneLists', '%s_symbol.gz') % (self.species), dtype=str, delimiter='\t').flatten()[1:].tolist()

        self.Convert(ensemblList, 'ensembl', 'hugo')
        self.Convert(symbolList, 'alias', 'hugo')

        return
    
    def Convert(self, genes, sourceType, targetType, onlineSearch = True, aggressiveSearch = False, returnUnknownString = True):
        
        class MyTypeError(Exception): 
            pass

        returnFlatFlag = False

        if (type(genes) is not list) and (type(genes) is not tuple):
            if (type(genes) is int) or (type(genes) is str):
                genes = [genes]
                returnFlatFlag = True
            else:
                raise MyTypeError('The only currently supported input types for "genes" are list, tuple, int and string')
        
        if onlineSearch:
            genesToFetch = set(genes).difference(self.conversionDict[sourceType][targetType].keys())

            if aggressiveSearch == False:
                genesToFetch = genesToFetch.difference(self.conversionDict[sourceType]['unknown'])

            if len(genesToFetch) > 0:
                self.Fetch(list(genesToFetch),sourceType)
        
        geneSet = set(self.conversionDict[sourceType][targetType].keys())
        
        if returnUnknownString:
            genes_converted = [ self.conversionDict[sourceType][targetType][gene] 
                                if gene in geneSet
                                else ('UNKNOWN (%s)' % gene) 
                                for gene in genes ]
        else:
            genes_converted = [ self.conversionDict[sourceType][targetType][gene] 
                                if gene in geneSet
                                else gene
                                for gene in genes ]
        
        if returnFlatFlag: 
            return genes_converted[0]
        else: 
            return genes_converted
    
    def Fetch(self, genes, sourceType):
        
        kwargs = {'species': self.species, 'fields': ['ensembl.gene', 'entrezgene', 'symbol', 'alias', 'retired']}

        if sourceType == 'entrez' or sourceType == 'retired':
            kwargs['scopes'] = ('entrezgene', 'retired')

            for gene in genes:
                assert type(gene) is int or type(gene) is long

        elif sourceType == 'ensembl':
            kwargs['scopes'] = 'ensembl.gene'

            for gene in genes:
                assert type(gene) is str
                assert gene[:4] == self.ensemblStr

        elif (sourceType == 'hugo') or (sourceType == 'alias'): 
            kwargs['scopes'] = ('symbol', 'alias')

            for gene in genes:
                assert type(gene) is str
                assert gene[:4] != self.ensemblStr
        
        genesSet = set(genes)
        queryList = mygene.MyGeneInfo().querymany(genes, **kwargs)
        
        #Loop through query results and fill in the dictionary
        for q,gene in zip(queryList,genes):
            #First check if gene was found by mygene
            if 'notfound' in q.keys():
                self.conversionDict[sourceType]['unknown'].add(gene)
            else:
                #Extract the different naming conventions from the query
                try:
                    entrez = q['entrezgene']

                    try:
                        retireds = q['retired']
                        if type(retireds) is int:
                            retireds = [retireds]
                        else:
                            retireds = list(retireds)
                        retireds.append(entrez)
                    except KeyError:
                        retireds = [entrez]
                except KeyError:
                    entrez = None
                    retireds = []

                try: 
                    hugo = q['symbol']
                except KeyError: 
                    hugo = None

                if str(q['query'])[:4] == self.ensemblStr:
                    ensembl = q['query']
                    ensembl_list = (ensembl,)
                else:
                    try:
                        ensembl = q['ensembl']

                        if type(ensembl) == dict:
                            ensembl = ensembl['gene']

                        if type(ensembl) == list:
                            ensembl_list = [list(i.values())[0] for i in ensembl]
                            ensembl = ensembl_list[0]

                        else:
                            ensembl_list = ()
                        
                    except KeyError:
                        ensembl = None
                        ensembl_list = ()
                try:
                    aliases = q['alias']

                    if type(aliases) != list: 
                        aliases = [aliases]

                    aliases.append(hugo)
                    aliases = tuple(set(aliases))

                except KeyError: 
                    aliases = ()
                
                #Map between all pairs of the three conventions (when possible)
                for source,sourceGene in zip(['entrez','hugo','ensembl'], [entrez,hugo,ensembl]):
                    for target,targetGene in zip(['entrez','hugo','ensembl'], [entrez,hugo,ensembl]):
                        if source is not None:
                            self.conversionDict[source]['known'].add(sourceGene)

                            if sourceGene in self.conversionDict[source]['unknown']:
                                self.conversionDict[source]['unknown'].remove(sourceGene)

                        if (source != target) and (sourceGene is not None) and (targetGene is not None):
                            self.conversionDict[source][target][sourceGene] = targetGene
                
                #Some searches identify many ensembl's for a given hugo/entrez query
                for ensembl_list_item in set(ensembl_list).difference(genesSet):
                    for target,targetGene in zip(['entrez','hugo'],[entrez,hugo]):
                        self.conversionDict['ensembl'][target][ensembl_list_item] = targetGene
                
                #Map all aliases individually to the other conventions
                for alias in aliases:
                    self.conversionDict['alias']['known'].add(alias)

                    if alias in self.conversionDict['alias']['unknown']:
                        self.conversionDict['alias']['unknown'].remove(alias)

                    for target,targetGene in zip(['entrez','hugo','ensembl'],[entrez,hugo,ensembl]):
                        if targetGene is not None:
                            self.conversionDict['alias'][target][alias] = targetGene
                
                #Map official hugo to a tuple of aliases
                if hugo is not None and len(aliases) > 0:
                    self.conversionDict['hugo']['alias'][hugo] = aliases
                
                #Map all retired entrez to official entrez and vice versa
                if entrez is not None:
                    self.conversionDict['entrez']['retired'][entrez] = retireds

                    for retired in retireds:
                        self.conversionDict['retired']['entrez'][retired] = entrez
                        self.conversionDict['retired']['known'].add(retired)

                        if retired in self.conversionDict['retired']['unknown']:
                            self.conversionDict['retired']['unknown'].remove(retired)
                
        #Fix alias dictionary to always assume that the input gene is the official
        #hugo.  This is important in the case of symbols like CRP, which is both an alias
        #of CSRP1 with EzID=1465 and an official hugo name with EzID=1401.  This final
        #loop guarantees
        #
        #   conversionDict['alias']['hugo']['CRP'] == 'CRP'
        #
        #Without this final loop, the alias dictionary would ERRONEOUSLY map the
        #official gene name CRP to CSRP1:
        #
        #   conversionDict['alias']['hugo']['CRP'] == 'CSRP1'
        
        for hugo in list(self.conversionDict['hugo']['entrez'].keys()) + \
                    list(self.conversionDict['hugo']['ensembl'].keys()) + \
                    list(self.conversionDict['hugo']['alias'].keys()):
            try: 
                entrez = self.conversionDict['hugo']['entrez'][hugo]
            except KeyError: 
                entrez = None

            try: 
                ensembl = self.conversionDict['hugo']['ensembl'][hugo]
            except KeyError: 
                ensembl = None

            for target, targetGene in zip(['entrez','hugo','ensembl'], [entrez,hugo,ensembl]):
                if targetGene is not None:
                    self.conversionDict['alias'][target][hugo] = targetGene
        
        write(self.conversionDict,self.dictDir)

        return






