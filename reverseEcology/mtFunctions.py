###############################################################################
# mtFunctions.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for working with metatranscriptomes.
################################################################################ 

# Import Python packages.
from Bio import SeqIO
import cobra
import networkx as nx
import os
import pandas as pd
import re

# These custom-written modules should have been included with the package
# distribution. 
import metadataFunctions as mf

################################################################################ 

# expressionMapping
# A function which identifies COGs and RPKM values associated with each
# metabolite in a network model. Inputs are the following directories:

# modelList: A list of models to perform this analysis on
# sampleList: A list of MT sets to perform this analysis on

# genomeModelDir: directory containing network and sbml models for individual
# genomes
# mergedModelDir: directory contaiing netowrk models and seed lists for 
# consensus genome models
# mergedDataDir: directory containing seed sets for each consensus genome
# genomeCogFolder: directory containing genes associated with each COG in a 
# genome
# rpkmFolder: directory containing RPKM values for each COG in a consensus
# genome
# taxonFile: File giving taxonomic information for individual genomes
# outFolder: where to put the output

# Output is a file MT-clade.out which gives RPKM values for each (metabolite, 
# COG) pairing

def expressionMapping(modelList, sampleList, genomeModelDir, mergedModelDir, 
                      mergedDataDir, genomeCogFolder, rpkmFolder, taxonFile,
                      outFolder):

    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
    
    lineageDict =  mf.importTaxonomy(taxonFile, 'Lineage')
    cladeDict =  mf.importTaxonomy(taxonFile, 'Clade')
    tribeDict =  mf.importTaxonomy(taxonFile, 'Tribe')
    sampleDict = dict(lineageDict.items() + cladeDict.items() + tribeDict.items())

    # Loop over each merged model...
    for sample in sampleList:
        
        print ('Processing sample: '+sample+': '+str(sampleList.index(sample)+1)+' of '+str(len(sampleList)))

        for curDir in modelList:
            
            print ('Processing genome: '+curDir+': '+str(modelList.index(curDir)+1)+' of '+str(len(modelList)))
            
            # Read its adjacency list representation
            mergedDiGraph = nx.read_adjlist(mergedModelDir+'/'+curDir+'/'+curDir+'RedAdjList.txt',
                                    create_using=nx.DiGraph())    
            
            # Create a dataFrame to store the results
            #    rpkmDF = pd.DataFrame(index=mergedDiGraph.nodes(), columns=['Seed', 'Cooperates with', 'COG', 'Direction', 'Avg RPKM', 'Norm RPKM'])
            summaryDF = pd.DataFrame(columns=['Metabolite', 'Seed', 'Cooperates with', 'COG', 'Direction', 'Avg RPKM', 'Norm RPKM'])
            
            # Retrieve the list of genomes associated with the merged model
            genomeList = sampleDict[curDir]
            
            # Read in the set of seed compounds
            seedDF = seedDF = pd.read_csv(mergedDataDir+'/'+curDir+'/'+curDir+'SeedCompounds.txt', 
                                          sep='\t', index_col=0, header=None, 
                                          names=['Metabolite', 'Metabolite', 'Weight'])
            
            # Loop over the metabolites in the model and create dictionaries for the 'Seed'
            #  and 'Cooperates' columns
            seedDict = dict.fromkeys(mergedDiGraph.nodes(), [])
            cooperatesDict = dict.fromkeys(mergedDiGraph.nodes(), [])
                
            for metab in mergedDiGraph.nodes():
                # Update 'Seed' column
                if metab in seedDF.index.tolist():
                    seedDict[metab] = 'Y'
                else:
                    seedDict[metab] = 'N'
                # Update 'Cooperates' column
                for innerDir in modelList:
                    if innerDir != curDir:
                        # Read in the dataframe of cooperative metabolites
                        coopSeedsDF = pd.read_csv(mergedDataDir+'/'+curDir+'/'+curDir+'-cooperates-'+innerDir+'.csv', sep=',', skiprows=1,
            usecols=[1, 2, 3],  index_col=0, header=None, names=['Index', 'Name', 'Inner Weight'])
                        if metab in coopSeedsDF.index:
                            if cooperatesDict[metab]:
                                cooperatesDict[metab] = cooperatesDict[metab]+','+innerDir
                            else:
                                cooperatesDict[metab] = innerDir
            for key in cooperatesDict:
                if not cooperatesDict[key]:
                    cooperatesDict[key] = "None"
            
            # Next, loop over genomes to create a link between the compounds in that 
            # genome and their COGs. This section will be run twice, once for reactions
            # which synthesize the compound, and once for reactions that consume the
            # compound.
            
            compoundToCogDictSynthesis = dict.fromkeys(mergedDiGraph.nodes(), [])
            compoundToCogDictConsumption = dict.fromkeys(mergedDiGraph.nodes(), [])
            
            for genome in genomeList:
            
                genomeEdgeMapping = pd.read_csv(genomeModelDir+'/'+genome+'/'+genome+'RxnEdges.txt', sep='\t', header=None, names=['Reaction'], index_col=[0,1])
                
                # Read in the SBML model containing GPRs
                genomeModel = cobra.io.read_sbml_model(genomeModelDir+'/'+genome+'/'+genome+'.xml')
            
                # Read in the COG dictionary
                cogDict = {}
                inFile = open(genomeCogFolder+'/'+genome+'COGs.txt', 'r')
                for line in inFile.readlines():
                    key = line.strip().split(',')[0]
                    # Reformat the key for compatability with COBRA
                    # Replace .genome. with .
                    key = re.sub('\.genome\.', '.', key)
                    # Replace . with underscore
                    key = re.sub('\.','_', key)
                    value = line.strip().split(',')[1]
                    cogDict[key] = value
                inFile.close()
                
                # For each metabolite in the consensus genome, find all sets of 
                # outward-pointing arcs - These are the "consumption" reactions
                # Also find all sets of inward-pointing arcs - these are the "synthesis"
                # reactions.
                # For each such, look up all associated reactions in the genome's adjacency list
                # For all reactions, look up all the associated genes in the model
                # For all genes, look up the proper COG
                for metab in mergedDiGraph.nodes():
                    # Look up all sink nodes
                    synthesizeNodeList = mergedDiGraph.predecessors(metab)
                    consumeNodeList = mergedDiGraph.successors(metab)
            
                    # Generate the list of associated reactions
                    synthesizeRxnList = []
                    consumeRxnList = []
            
                    for node in synthesizeNodeList:
                        if (node, metab) in genomeEdgeMapping.index.tolist():
                            synthesizeRxnList.append(genomeEdgeMapping.loc[node, metab]['Reaction'].tolist())
                    # Reduce to unique entries
                    synthesizeRxnList = [rxn for rxnList in synthesizeRxnList for rxn in rxnList]
                    synthesizeRxnList = list(set(synthesizeRxnList))
            
                    for node in consumeNodeList:
                        if (metab, node) in genomeEdgeMapping.index.tolist():
                            consumeRxnList.append(genomeEdgeMapping.loc[metab, node]['Reaction'].tolist())
                    # Reduce to unique entries
                    consumeRxnList = [rxn for rxnList in consumeRxnList for rxn in rxnList]
                    consumeRxnList = list(set(consumeRxnList))
                            
                    # Generate the list of associated genes
                    synthesizeGeneList = []
                    for reaction in synthesizeRxnList:
                        for gene in genomeModel.reactions.get_by_id(reaction).genes:
                            if gene.id != 'Unknown':                    
                                synthesizeGeneList.append(gene.id)
                    # Reduce to unique entries
                    synthesizeGeneList = list(set(synthesizeGeneList))
                
                    consumeGeneList = []
                    for reaction in consumeRxnList:
                        for gene in genomeModel.reactions.get_by_id(reaction).genes:
                            if gene.id != 'Unknown':                    
                                consumeGeneList.append(gene.id)
                    # Reduce to unique entries
                    consumeGeneList = list(set(consumeGeneList))
                    
                    # Generate the list of associated COGs
                    synthesizeCogList = []
                    for gene in synthesizeGeneList:
                        synthesizeCogList.append(cogDict[gene])
                        # Reduce to unique entries
                    synthesizeCogList = list(set(synthesizeCogList))
                    
                    consumeCogList = []
                    for gene in consumeGeneList:
                        consumeCogList.append(cogDict[gene])
                        # Reduce to unique entries
                    consumeCogList = list(set(consumeCogList))
                
                    # Update compoundToCogDict with new COG elements
                    compoundToCogDictSynthesis[metab] = compoundToCogDictSynthesis[metab] + synthesizeCogList
                    compoundToCogDictConsumption[metab] = compoundToCogDictConsumption[metab] + consumeCogList
            
            # Now that the dict has been fully populated, reduce keys to lists of unique values
            for metab in compoundToCogDictSynthesis.keys():
                compoundToCogDictSynthesis[metab] = list(set(compoundToCogDictSynthesis[metab]))                
            for metab in compoundToCogDictConsumption.keys():
                compoundToCogDictConsumption[metab] = list(set(compoundToCogDictConsumption[metab]))             
            
            # Read in the RPKM data
            rpkmDF = pd.read_csv(rpkmFolder+'/'+sample+'.norm', sep=',', index_col=0)
            
            # Now we need to populate the dataframe
            for metab in mergedDiGraph.nodes():
                synthesizeCogList = compoundToCogDictSynthesis[metab]
                for cog in synthesizeCogList:
                    summaryDF = summaryDF.append({'Metabolite': metab, 'Seed': seedDict[metab],
                                                  'Cooperates with': cooperatesDict[metab],
                                                  'COG': cog, 'Direction': 'Synthesis',
                                                  'Avg RPKM': rpkmDF.loc[curDir+'-'+cog, 'Avg RPKM'],
                                                  'Norm RPKM': rpkmDF.loc[curDir+'-'+cog, 'Norm RPKM'],
                                                  }, ignore_index=True)
                consumeCogList = compoundToCogDictConsumption[metab]
                for cog in consumeCogList:
                    summaryDF = summaryDF.append({'Metabolite': metab, 'Seed': seedDict[metab],
                                                  'Cooperates with': cooperatesDict[metab],
                                                  'COG': cog, 'Direction': 'Consumption',
                                                  'Avg RPKM': rpkmDF.loc[curDir+'-'+cog, 'Avg RPKM'],
                                                  'Norm RPKM': rpkmDF.loc[curDir+'-'+cog, 'Norm RPKM'],
                                                  }, ignore_index=True)
            
            summaryDF.to_csv(outFolder+'/'+sample+'-'+curDir+'.out', index=False)
            
################################################################################ 

# seedsToGenes
# A function which identifies genes associated with each seed metabolite in a 
# network model. Inputs are the following directories:

# modelList: A list of models to perform this analysis on

# modelDir: directory containing network and sbml models for individual
# genomes
# dataDir: directory containing seed sets for each consensus genome
# genbankFolder: folder containing genbank files of individual genomes

# Output is a file genome.out which lists the reaction and annotation for each
# (metabolite, gene) pairing

def seedsToGenes(modelList, modelDir, dataDir, genbankFolder):
        
    # Loop over each genome...
    i = 1
    for curDir in modelList:
        print "Processing model "+str(i)+" of "+str(len(modelList))
        # Read its adjacency list representation
        mergedDiGraph = nx.read_adjlist(modelDir+'/'+curDir+'/'+curDir+'RedAdjList.txt',
                                create_using=nx.DiGraph())    
        
        # Create a dataFrame to store the results
        #    rpkmDF = pd.DataFrame(index=mergedDiGraph.nodes(), columns=['Seed', 'Cooperates with', 'COG', 'Direction'])
        summaryDF = pd.DataFrame(columns=['Metabolite', 'Seed', 'Cooperates with', 'Gene', 'Reaction', 'Annotation', 'Direction'])
        # Read in the set of seed compounds
        seedDF = seedDF = pd.read_csv(dataDir+'/'+curDir+'/'+curDir+'SeedCompounds.txt', 
                                      sep='\t', index_col=0, header=None, 
                                      names=['Metabolite', 'Metabolite', 'Weight'])
        
        # Loop over the metabolites in the model and create dictionaries for the 'Seed'
        #  and 'Cooperates' columns
        seedDict = dict.fromkeys(mergedDiGraph.nodes(), [])
        cooperatesDict = dict.fromkeys(mergedDiGraph.nodes(), [])
            
        for metab in mergedDiGraph.nodes():
            # Update 'Seed' column
            if metab in seedDF.index.tolist():
                seedDict[metab] = 'Y'
            else:
                seedDict[metab] = 'N'
            # Update 'Cooperates' column
            for innerDir in modelList:
                if innerDir != curDir:
                    # Read in the dataframe of cooperative metabolites
                    coopSeedsDF = pd.read_csv(dataDir+'/'+curDir+'/'+curDir+'-cooperates-'+innerDir+'.csv', sep=',', skiprows=1,
        usecols=[1, 2, 3],  index_col=0, header=None, names=['Index', 'Name', 'Inner Weight'])
                    if metab in coopSeedsDF.index:
                        if cooperatesDict[metab]:
                            cooperatesDict[metab] = cooperatesDict[metab]+','+innerDir
                        else:
                            cooperatesDict[metab] = innerDir
        for key in cooperatesDict:
            if not cooperatesDict[key]:
                cooperatesDict[key] = "None"
        
        # Next, loop over genomes to create a link between the compounds in that 
        # genome and their COGs. This section will be run twice, once for reactions
        # which synthesize the compound, and once for reactions that consume the
        # compound.
        
        compoundToGeneDictSynthesis = dict.fromkeys(mergedDiGraph.nodes(), [])
        compoundToGeneDictConsumption = dict.fromkeys(mergedDiGraph.nodes(), [])
        
        genomeEdgeMapping = pd.read_csv(modelDir+'/'+curDir+'/'+curDir+'RxnEdges.txt', sep='\t', header=None, names=['Reaction'], index_col=[0,1])
        
        # Read in the SBML model containing GPRs
        genomeModel = cobra.io.read_sbml_model(modelDir+'/'+curDir+'/'+curDir+'.xml')
        
        # For each metabolite in the consensus genome, find all sets of 
        # outward-pointing arcs - These are the "consumption" reactions
        # Also find all sets of inward-pointing arcs - these are the "synthesis"
        # reactions.
        # For each such, look up all associated reactions in the genome's adjacency list
        # For all reactions, look up all the associated genes in the model
        # For all genes, look up the proper COG
        for metab in mergedDiGraph.nodes():
            # Look up all sink nodes
            synthesizeNodeList = mergedDiGraph.predecessors(metab)
            consumeNodeList = mergedDiGraph.successors(metab)
    
            # Generate the list of associated reactions
            synthesizeRxnList = []
            consumeRxnList = []
    
            for node in synthesizeNodeList:
                if (node, metab) in genomeEdgeMapping.index.tolist():
                    synthesizeRxnList.append(genomeEdgeMapping.loc[node, metab]['Reaction'].tolist())
            # Reduce to unique entries
            synthesizeRxnList = [rxn for rxnList in synthesizeRxnList for rxn in rxnList]
            synthesizeRxnList = list(set(synthesizeRxnList))
    
            for node in consumeNodeList:
                if (metab, node) in genomeEdgeMapping.index.tolist():
                    consumeRxnList.append(genomeEdgeMapping.loc[metab, node]['Reaction'].tolist())
            # Reduce to unique entries
            consumeRxnList = [rxn for rxnList in consumeRxnList for rxn in rxnList]
            consumeRxnList = list(set(consumeRxnList))
                    
            # Generate the list of associated genes
            synthesizeGeneRxnList = []
            for reaction in synthesizeRxnList:
                for gene in genomeModel.reactions.get_by_id(reaction).genes:
                    if gene.id != 'Unknown':                    
                        synthesizeGeneRxnList.append([gene.id.replace('_CDS_', '.CDS.'), reaction])
    
            consumeGeneRxnList = []
            for reaction in consumeRxnList:
                for gene in genomeModel.reactions.get_by_id(reaction).genes:
                    if gene.id != 'Unknown':                    
                        consumeGeneRxnList.append([gene.id.replace('_CDS_', '.CDS.'), reaction])
        
            # Update compoundToCogDict with new COG elements
            compoundToGeneDictSynthesis[metab] = compoundToGeneDictSynthesis[metab] + synthesizeGeneRxnList
            compoundToGeneDictConsumption[metab] = compoundToGeneDictConsumption[metab] + consumeGeneRxnList
        
        # Establish a dataframe linking gene IDs to functions. This segment reads in
        # the genbank file and creates a dict linking CDS to gene IDs for each gene.
        geneAnnotDict = {}
        for contig in SeqIO.parse(open(genbankFolder+'/'+curDir+'.gbk', "r"), 'genbank'):
            for feature in contig.features:
                if feature.type == 'CDS':
                    gene = feature.qualifiers['gene'][0].replace('.genome.', '.')
                    annot = feature.qualifiers['function'][0]
                    geneAnnotDict[gene] = annot
    
        # Now we need to populate the dataframe
        for metab in mergedDiGraph.nodes():
            synthesizeGeneList = compoundToGeneDictSynthesis[metab]
            for geneRxnPair in synthesizeGeneList:
                summaryDF = summaryDF.append({'Metabolite': metab, 'Seed': seedDict[metab],
                                              'Cooperates with': cooperatesDict[metab],
                                              'Gene': geneRxnPair[0], 'Reaction': geneRxnPair[1],
                                              'Annotation':  geneAnnotDict[geneRxnPair[0]],                                          
                                              'Direction': 'Synthesis',
                                              }, ignore_index=True)
            consumeGeneList = compoundToGeneDictConsumption[metab]
            for geneRxnPair in consumeGeneList:
                summaryDF = summaryDF.append({'Metabolite': metab, 'Seed': seedDict[metab],
                                              'Cooperates with': cooperatesDict[metab],
                                              'Gene': geneRxnPair[0], 'Reaction': geneRxnPair[1],
                                              'Annotation':  geneAnnotDict[geneRxnPair[0]],                                          
                                              'Direction': 'Consumption',
                                              }, ignore_index=True)
        
        i = i + 1
        
        summaryDF.to_csv(modelDir+'/'+curDir+'/'+curDir+'.out', index=False)     
    return
    
################################################################################ 

# mergedRedAdjListToCogList
# A function which identifies the COGs associated with the reduced adjacency
# list for a clade. These represent the COGs from which seed sets are computed;
# the remaining COGs should be investigated for interesting activities which
# did not get captured by the seed set analysis.

def mergedRedAdjListToCogList(cladeList, genomeModelDir, genomeCogFolder, taxonFile, mergedModelDir):

    ################################################################################
    ### For each clade, find the list of COGs associated with the reduced ajacency
    ### list from which seed compounds were computed. This requires a complicated
    ### mapping:
    ###   Read in the list of (cmpd, cmpd) arcs which make up the redAdjList
    ###   For each genome in the clade, read its RxnEdges mapping, to link these
    ###   arcs to reactions
    ###   For each genome in the clade, read its COBRA model, to link these
    ###   reactions to genes
    ###   For each genome in the clade, read its genomeCOGs file, to link these
    ###   genes to COGs
    ### This will give three dictionaries, linking: arcs to reactions, reactions
    ### to genes, and genes to COGs.
    ### Then, iterate once more over the redAdjList, aggregating the COGs
    ################################################################################
    
    # Create a dict of genomes for each clade
    cladeDict = mf.importTaxonomy(taxonFile, 'Clade')
    
    # Loop over each clade to perform the analysis described above
    for clade in cladeList:
    
        # Create empty dictionaries to store links between arcs, reactions, genes,
        # and COGs
        arcToRxnDict = {}
        rxnToGeneDict = {}
        geneToCogDict = {}
        
        for genome in cladeDict[clade]:
            # Create the arcToRxnDict
            # Because an arc may be associated with multiple reactions, we need
            # a check
            inFile = open(genomeModelDir+'/'+genome+'/'+genome+'RxnEdges.txt', 'r')
            for line in inFile.readlines():
                lineArray = line.strip().split('\t')
                # Check to see if the arc has been detected before. If yes, append
                # the reaction. Otherwise, create a new dictionary entry.
                arc = lineArray[0]+'-'+lineArray[1]
                if arc in arcToRxnDict.keys():
                    arcToRxnDict[arc].append(lineArray[2])
                else:
                    arcToRxnDict[arc] = [lineArray[2]]
     
            # Create the rxnToGeneDict
            # Because a reaction may be associated with multiple genes, we need
            # a check
            genomeModel = cobra.io.read_sbml_model(genomeModelDir+'/'+genome+'/'+genome+'.xml')
            for curRxn in genomeModel.reactions:
                for gene in curRxn.genes:
                    if gene.id != 'Unknown':
                        # Check to see if the reaction has been detected before. If yes, 
                        # append the gene. Otherwise, create a new dictionary entry.
                        if curRxn.id in rxnToGeneDict.keys():
                            rxnToGeneDict[curRxn.id].append(gene.id)
                        else:
                            rxnToGeneDict[curRxn.id] = [gene.id]
                    
            # Create the geneToCogDict
            inFile = open(genomeCogFolder+'/'+genome+'COGs.txt', 'r')
            for line in inFile.readlines():
                key = line.strip().split(',')[0]
                # Reformat the key for compatability with COBRA
                # Replace .genome. with .
                key = re.sub('\.genome\.', '.', key)
                # Replace . with underscore
                key = re.sub('\.','_', key)
                value = line.strip().split(',')[1]
                geneToCogDict[key] = value
            inFile.close()
            
        # Reduce dictionary lists to their unique elements
        for key in arcToRxnDict.keys():
            arcToRxnDict[key] = list(set(arcToRxnDict[key]))
    
        # Read in the redAdjList as a diGraph
        mergedDiGraph = nx.read_adjlist(mergedModelDir+'/'+clade+'/'+clade+'RedAdjList.txt', create_using=nx.DiGraph())
        
        arcToCogDict = {}
    
        for arc in mergedDiGraph.edges():
            # Convert arc to representation used above
            arc = list(arc)
            arc = arc[0]+'-'+arc[1]
            
            # Iterate through the appropriate dicts to retrieve the COGs for each
            # arc
            rxnList = arcToRxnDict[arc]
            for rxn in rxnList:
                geneList = rxnToGeneDict[rxn]
                for gene in geneList:
                    cog = geneToCogDict[gene]
                    if arc in arcToCogDict.keys():
                        arcToCogDict[arc].append(cog)
                    else:
                        arcToCogDict[arc] = [cog]
                        
        # Reduce dictionary lists to their unique elements
        for key in arcToCogDict.keys():
            arcToCogDict[key] = list(set(arcToCogDict[key]))
            
        # Loop over dictionary to generate a list of all COGs
        cogList = []
        for key in arcToCogDict.keys():
            cogList.append(arcToCogDict[key])
        
        cogList = [item for sublist in cogList for item in sublist]
    
        # Reduce the cogList to its unique entries
        cogList=list(set(cogList))
            
        # Write the results to file
        outFile = open(mergedModelDir+'/'+clade+'/'+clade+'COGs.txt', 'w')
        for cog in cogList:
            outFile.write(cog+'\n')
        outFile.close()
    
    return