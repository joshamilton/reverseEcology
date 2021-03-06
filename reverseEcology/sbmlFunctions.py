###############################################################################
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for manipulating SBML files
################################################################################

# Import python modules
import cobra
import cobra.core.Formula
import collections
import fileinput
import numpy as np
import os
import pandas as pd
import re

# Import custom Python modules 
import metadataFunctions as mf

# Define path for data included in the package
dataPath = os.path.dirname(os.path.abspath(__file__))+'/packageData'

################################################################################
    
# dirListToAdjacencyList
# This function iterates over a list of genome directories and converts each
# genome scale model from an SBML file to an adjacency list. Adjacency lists 
# for each genome-scale model are written as text files in each genome 
# directory. Summary statistics about each graph are written in the
# summaryStatsDir as well.

def dirListToAdjacencyList(dirList, processedDataDir, summaryStatsDir):

    numSubDir = len(dirList)

# Create an array to store summary statistics. The array has three integer
# columns, which will contain the number of genes, metabolites, and reactions
# in the SBML file. 
    modelStatArray = np.empty([numSubDir, 3], dtype = int)

# Create a file to record the summary statistics.
    modelFile = open(summaryStatsDir+'/'+'ModelStatistics.txt', 'w')
    modelFile.write('Model,Genes,Metabolites,Reactions\n')

# Iterate over the list of genome directories. For each genome, read in the
# SBML file and update the 'description' field with the genome name. The number
# of genes, metabolites, and reactions in the SBML file is recorded in the
# 'modelStatArray' and written to 'modelFile.' Finally, the genome-scale model
# is converted to an adjacency list and written to file.
    count = 0
    print 'Converting SBML files to adjacency lists'

# Create an empty dictionary to store metabolite IDs and names
    namesDict = {}
    
    for curDir in dirList:
# Read in SBML file    
        model = cobra.io.read_sbml_model(processedDataDir+'/'+curDir+'/'+curDir+'.xml')

# Create dictionary of metabolite names
        for metab in model.metabolites:
            namesDict[metab.id] = metab.name

# Update description field
        model.id = curDir;

# Read model statistics by invoking sbmlFunctions.getModelStats
        modelStatArray[count:] = getModelStats(model)
        modelFile.write('%s,%i,%i,%i\n' % (processedDataDir+'/'+curDir, modelStatArray[count,0], 
                                    modelStatArray[count,1], 
                                    modelStatArray[count, 2] ) )

# Create adjacency list and write to file
        adjacencyListFromModel(model, processedDataDir)
        reactionEdgesFromModel(model, processedDataDir)
        count = count + 1

# Close files containing summary data
    modelFile.close()

    return modelStatArray

################################################################################

# adjacencyListFromModel
# Function to convert a cobrapy model object to an adjaceny list. Also writes
# adjacency list to file. For each reaction in the model, creates an edge
# between all (reactant, product) pairs. If a reaction is reversible, also
# creates edges between all (product, reactant) pairs.
# Input: cobrapy model object, model directory
# Output: None.

def adjacencyListFromModel(model, processedDataDir):

# Establish a file for the adjacency list
    myFile = open(processedDataDir+'/'+model.id+'/'+model.id+'AdjList.txt', 'w')

# For each reaction, loop over the reactants. For each reactant, loop over the 
# reaction products and create an edge between the reactant and products. If a 
# reaction is reversible, repeat the process in reverse, creating an edge
# between each product and reactant.
    for myRxn in model.reactions:
        for myReactant in myRxn.reactants:
            myFile.write(myReactant.id+'\t')
            for myProduct in myRxn.products:
                myFile.write(myProduct.id+'\t')
            myFile.write('\n')
        if myRxn.reversibility == True:
            for myProduct in myRxn.products:
                myFile.write(myProduct.id+'\t')
                for myReactant in myRxn.reactants:
                    myFile.write(myReactant.id+'\t')
                myFile.write('\n')
    myFile.close()
    return
    
################################################################################

# reactionEdgesFromModel
# Function to convert a cobrapy model object to a list of (source, sink) pairs
# for reaction in the model. For each reaciton, creates an edge between all 
# (reactant, product) pairs and indicates the reaction. If a reaction is 
# reversible, also creates edges between all (product, reactant) pairs.
# Input: cobrapy model object, model directory
# Output: None.

def reactionEdgesFromModel(model, processedDataDir):

# Establish a file for the adjacency list
    myFile = open(processedDataDir+'/'+model.id+'/'+model.id+'RxnEdges.txt', 'w')

# For each reaction, loop over the reactants. For each reactant, loop over the 
# reaction products and create an edge between the reactant and products. If a 
# reaction is reversible, repeat the process in reverse, creating an edge
# between each product and reactant. Also record reaction associated with each
# edge.
    for myRxn in model.reactions:
        for myReactant in myRxn.reactants:
            for myProduct in myRxn.products:
                myFile.write(myReactant.id+'\t')
                myFile.write(myProduct.id+'\t')
                myFile.write(myRxn.id+'\n')
        if myRxn.reversibility == True:
            for myProduct in myRxn.products:
                for myReactant in myRxn.reactants:
                    myFile.write(myProduct.id+'\t')
                    myFile.write(myReactant.id+'\t')
                    myFile.write(myRxn.id+'\n')
    myFile.close()
    return
    
################################################################################

# getModelStats
# Function to retrieve statistics about a model
# Input: cobrapy model object of a genome-scale model
# Output: array with three integer columns, containing the number of genes, 
# metabolites, and reactions in the model

def getModelStats(model):
    statRow = [0]*3
    statRow[0] = len(model.genes)
    statRow[1] = len(model.metabolites)
    statRow[2] = len(model.reactions)
    return statRow

################################################################################

# Draft reconstructions from Kbase require some post-processing. This script 
# does several important things:
# 1. Reformat gene locus tags
# 2. Remove biomass, exchange, spontaneous, DNA/RNA biosynthesis reactions and 
# their corresponding genes
# 3. Import metabolite formulas
# 4. Check mass- and charge-balancing of reactions in the reconstruction
# 5. Remove trailing 0s from reaction and metabolite names

# The post-processing has a major shortcoming. When KBase detects that one or 
# more subunits of a complex are present, it creates a "full" GPR by adding 
# 'Unknown' genes for the other subunits. CobraPy currently lacks functions to 
# remove the genes. As such, these model should not be used to perform any 
# simulations which rely on GPRs.

# Each model should be in its own directory in the 'RawModelFiles' folder. Both
# SBMl and TSV versions from KBase are required.

# As output, the code returns processed SBML files in the 'processedDataDir'
# folder. Also returns a summary of the model sizes, in the 'summaryStatsDir'
# folder.

def processSBMLforRE(rawModelDir, processedDataDir, summaryStatsDir):

    # # Check that folders exist and create them if necessary
    if not os.path.exists(processedDataDir):
        os.makedirs(processedDataDir)
    if not os.path.exists(summaryStatsDir):
        os.makedirs(summaryStatsDir)
    
    # Import the list of models
    dirList = mf.getDirList(rawModelDir)
    numSubDir = len(dirList)
    
    # Import the list of metabolies to revise
    metabFormDict = {}
    metabChargeDict = {}
    
    with open(dataPath+'/newFormulaDict.txt') as inFile:
        for line in inFile:
           (key, value) = line.split('\t')
           metabFormDict[key] = value
           
    with open(dataPath+'/newChargeDict.txt') as inFile:
        for line in inFile:
           (key, value) = line.split('\t')
           metabChargeDict[key] = value
           
    # Create an array to store results
    # Columns: genes, metabs, rxns, balanced (binary)
    modelSizeDF = pd.DataFrame(index = dirList, columns=['Genes', 'Metabolites', 'Reactions', 'Balanced'])
    
    # Intialize a counter
    count = 1
    for curDir in dirList:

    # Print the subdirectory name
        print 'Processing model ' + curDir + ', ' + str(count) + ' of ' + str(numSubDir)
    
    # Import metabolite charges
        cpdData = pd.read_csv(rawModelDir+'/'+curDir+'/'+curDir+'Compounds.tsv', delimiter='\t', index_col=0)
    
    ################################################################################                   
    
    # Before reading in the SBML file, update gene loci so they read in properly
    # In the old models ...
    # KBase generates gene loci of the form kb|g.######.CDS.###
    # Transform them to be of the form curDir_CDS_###
    
    # In the new models ...
    # KBase generates gene loci of the form curDir.genome.CDS.###
    # Transform them to be of the form curDir_CDS_###
    
        for myLine in fileinput.FileInput(rawModelDir+'/'+curDir+'/'+curDir+'.xml', inplace = True):
    #            print re.sub('kb\|g\.\d+\.CDS\.(\d+)', curDir+'_CDS_\g<1>', myLine).strip()
            print re.sub(curDir+'\.genome\.CDS\.(\d+)', curDir+'_CDS_\g<1>', myLine).strip()
        fileinput.close()
    
    ################################################################################                   
    
    # Read in model from SBML
        model = cobra.io.read_sbml_model(rawModelDir+'/'+curDir+'/'+curDir+'.xml')
    
    ################################################################################                   
    
    # Remove undesired reactions, including:
        badRxnList = []
    #        print 'Removing bad reactions'
        for curRxn in model.reactions:
        # Exchange reactions
            if re.match('EX_', curRxn.id):
                badRxnList.append(curRxn)
        # Reactions w/o GPRs
            elif len(curRxn.gene_reaction_rule) == 0:
                badRxnList.append(curRxn)
        # Protein, DNA, and RNA synthesis
            elif curRxn.id == 'rxn13782_c0' or curRxn.id == 'rxn13783_c0' or curRxn.id == 'rxn13784_c0':
                badRxnList.append(curRxn)
        # Spontaneous reactions, whose GPR is fully 'unknown'
            elif curRxn.gene_reaction_rule == 'Unknown':
                badRxnList.append(curRxn)     
        # Transport reactions, based on keywords
            elif re.search('transport', curRxn.name) or re.search('permease', curRxn.name) or re.search('symport', curRxn.name) or re.search('diffusion', curRxn.name) or re.search('excretion', curRxn.name) or re.search('export', curRxn.name) or re.search('secretion', curRxn.name) or re.search('uptake', curRxn.name) or re.search('antiport', curRxn.name):
                badRxnList.append(curRxn)
        # Transport reactions which don't get picked up based on keywords
            elif curRxn.id == 'rxn05226_c0' or curRxn.id == 'rxn05292_c0' or curRxn.id == 'rxn05305_c0' or curRxn.id == 'rxn05312_c0' or curRxn.id == 'rxn05315_c0' or curRxn.id == 'rxn10945_c0' or curRxn.id == 'rxn10116_c0':
                badRxnList.append(curRxn)
        # Transport reactions, where a metabolite with the same ID is on both sides
            metabList = []
            for curMetab in curRxn.metabolites:
                metabList.append(re.sub('_[a-z]\d', '', curMetab.id))
            # Count number of appearences each list element
            metabList = [metab for metab in metabList if metab != 'cpd00067']
            if len(metabList) > 0:        
                metabCounter = collections.Counter(metabList)
                # Sort by number of appearances, find appearances of most common
                # Any value > 1 means the metablite is on both sides
                if metabCounter.most_common()[0][1] > 1:
                    badRxnList.append(curRxn)
    
        badRxnList = list(set(badRxnList))
        model.remove_reactions(badRxnList, delete=True, remove_orphans=True)                        
    
        print 'The remaining extracellular metabolites are:'
        for curMetab in model.metabolites:
            if re.search('_e0', curMetab.id):
                print curMetab.id
    
    ################################################################################                   
    
    # Update the metabolite formulas
    #        print 'Updating metabolite formulas'
        for curMetab in model.metabolites:
        # Retrieve the metabolite name w/o compartment info
        # Look up the appropriate index in the cpdData DF
            curMetab.formula = cobra.core.Formula.Formula(cpdData.loc[re.sub('_[a-z]\d', '', curMetab.id)][1])
            curMetab.formula.id = cpdData.loc[re.sub('_[a-z]\d', '', curMetab.id)][1]
    
    ################################################################################                   
    
    # Check mass- and charge- balancing   
        imbalCounter = 0
    #        print 'Correcting mass- and charge-balancing'
    
    # Check for reactions known to be imbalanced and manually correct them
    # If metabolite cpd03422 exists, update its charge to +1
        for curMetab in model.metabolites:
            if curMetab.id in metabFormDict.keys():
                curMetab.formula = cobra.core.Formula.Formula(metabFormDict[curMetab.id])
                curMetab.formula.id = metabFormDict[curMetab.id]
                curMetab.charge = int(metabChargeDict[curMetab.id])
    
        for curRxn in model.reactions:
    # If reaction rxn05893 exists, update its stoichiometry
            if curRxn.id == 'rxn05893_c0':
                curRxn.reaction = '4.0 cpd00001_c0 + 2.0 cpd00013_c0 + 6.0 cpd11621_c0 <=> 16.0 cpd00067_c0 + 2.0 cpd00075_c0 + 6.0 cpd11620_c0'         
    # If reaction rxn07295 exists, update its stoichiometry
            if curRxn.id == 'rxn07295_c0':
                curRxn.reaction = 'cpd00007_c0 + cpd00033_c0 <=> cpd00025_c0 + 3.0 cpd00067_c0 + cpd14545_c0'
    # If reaction rxn08808 exists, update its stoichiometry        
            elif curRxn.id == 'rxn08808_c0':
                curRxn.reaction = 'cpd00001_c0 + cpd15341_c0 <=> cpd00067_c0 + cpd00908_c0 + cpd01080_c0'
    # If reaction rxn12822 exists, update its stoichiometry        
            elif curRxn.id == 'rxn12822_c0':
                curRxn.reaction = '2.0 cpd00023_c0 + cpd11621_c0 <=> cpd00024_c0 + cpd00053_c0 + 2.0 cpd00067_c0 + cpd11620_c0'
    # If reaction rxn12822 exists, update its stoichiometry        
    
    # Heuristic for proton balancing
        for curRxn in model.reactions:
            imbalDict = curRxn.check_mass_balance()
            if len(imbalDict) != 0:
            # If imbalancing due to protons alone, correct it
                if 'H' in imbalDict and 'charge' in imbalDict and imbalDict['H'] == imbalDict['charge']:
                    curRxn.add_metabolites({model.metabolites.get_by_id('cpd00067_c0'): -1*imbalDict['H']})
                else:
                    imbalCounter = imbalCounter + 1                
                    print 'Reaction ' + str(curRxn.id) + ' remains unbalanced'      
                
    # Inform of results
        if imbalCounter != 0:
            modelSizeDF.loc[curDir][3] = 0
        else:
            modelSizeDF.loc[curDir][3] = 1
            print 'All reactions are balanced'              
        
    ################################################################################            
    
    ### Update names to remove trailing zeros
        for curComp in model.compartments:
            model.compartments[curComp] = re.sub('_\d', '', model.compartments[curComp])
            model.compartments[re.sub('\d', '', curComp)] = model.compartments.pop(curComp)
        
        for curMetab in model.metabolites:
            curMetab.id = re.sub('\d$', '', curMetab.id)
            curMetab.name = re.sub('_[a-z]\d$', '', curMetab.name)
            curMetab.compartment = re.sub('\d$', '', curMetab.compartment)
        
        for curRxn in model.reactions:
            curRxn.id = re.sub('\d$', '', curRxn.id)
            curRxn.name = re.sub('_[a-z]\d$', '', curRxn.name)
                
    ################################################################################                   
    
    # Store the model properties in the array, write the model output, and increase the counter
        modelSizeDF.loc[curDir][0] = len(model.genes)
        modelSizeDF.loc[curDir][1] = len(model.metabolites)
        modelSizeDF.loc[curDir][2] = len(model.reactions)
    
    # Perform final write to file
    # Check that output dir exists and create if necessary    
        if not os.path.exists(processedDataDir+'/'+curDir):
            os.makedirs(processedDataDir+'/'+curDir)
    
    #        print 'Writing to file'
        cobra.io.write_sbml_model(model, processedDataDir+'/'+curDir+'/'+curDir+'.xml')
        count = count + 1
    
    # Write the results to file
    modelSizeDF.to_csv(summaryStatsDir+'/modelStats.tsv', sep='\t')    

    return
    
################################################################################

# Prior to reverse ecology analysis, we "prune" the network topology to make
# the arcs in the directed graph more "physiologically realistic." The criteria
# we use are outlined in 

# Ma, H., & Zeng, A. P. (2003). Reconstruction of metabolic networks from 
# genome data and analysis of their global structure for various organisms. 
# Bioinformatics, 19(2) 270-277.

# Briefly, we first remove all pairs of currency metabolites involved in
# proton or functional group transfer. Then, we remove additional singletom
# metabolites (protons and the like.)

def pruneCurrencyMetabs(modelDir, summaryStatsDir):
    
    # Import the list of models
    dirList = mf.getDirList(modelDir)
    
    # Create an array to store results
    # Columns: genes, metabs, rxns, balanced (binary)
    modelSizeDF = pd.DataFrame(index = dirList, columns=['Genes', 'Metabolites', 'Reactions'])
    
    # Intialize a counter
    count = 1
    
    # Process each model...
    for curDir in dirList:

    # Read in model from SBML
        model = cobra.io.read_sbml_model(modelDir+'/'+curDir+'/'+curDir+'.xml')
    
    # Write the original model to a text file for inspection.
    #    with open(modelDir+'/'+curDir+'/'+curDir+'Original.txt', 'w') as outFile:
    #        for curRxn in model.reactions:
    #            outFile.write(curRxn.id+'\t'+curRxn.build_reaction_string(use_metabolite_names=True)+'\n')
    
    #############################
                    
    # Read in the list of bad metabolite pairs
        pairMetabList = []
        with open(dataPath+'/currencyRemovePairs.txt') as myFile:
            for line in myFile:
                pairMetabList.append(line.strip().split('\t'))
    
    # Remove all pairs of metabolites
        for curRxn in model.reactions:
            for pair in pairMetabList:
                # List metabolites in the reaction. Needs to be inside loop b/c reaction metabolites will change
                metabList = []
                for curMetab in curRxn.metabolites:
                    metabList.append(curMetab.id)
                # If both members of the pair participate in the reaction, drop both from the reaction
                if set(pair) <= set(metabList):
                    for metab in pair:
                        curRxn.pop(model.metabolites.get_by_id(metab))
    
    #############################
    
    # Read in the list of aminotransfer metabolite pairs
        pairAminoList = []
        with open(dataPath+'/currencyAminoPairs.txt') as myFile:
            for line in myFile:
                pairAminoList.append(line.strip().split('\t'))
    
    # Remove all pairs of metabolites
        for curRxn in model.reactions:
            for pair in pairAminoList:
                # List metabolites in the reaction. Needs to be inside loop b/c reaction metabolites will change
                metabList = []
                for curMetab in curRxn.metabolites:
                    metabList.append(curMetab.id)
                # If both members of the pair participate in the reaction, drop both from the reaction
                if (set(pair) <= set(metabList)) and ('cpd00013_c' not in metabList):
                    for metab in pair:
                        curRxn.pop(model.metabolites.get_by_id(metab))
                    
    #############################
    
    # Read in the list of bad metabolite singletons
        with open(dataPath+'/currencyRemoveSingletons.txt') as myFile:
            singletonMetabList = myFile.read().splitlines()
    
    # Remove all bad metabolites
        for curMetab in singletonMetabList:
            if model.metabolites.has_id(curMetab):
                model.metabolites.get_by_id(curMetab).remove_from_model(method='subtractive')
    
    # Prune the model, dropping empty reactions and those with only a subtrate or product
        cobra.manipulation.delete.prune_unused_reactions(model)
        for curRxn in model.reactions:
            if len(curRxn.reactants) == 0 or len(curRxn.products) == 0:
                curRxn.remove_from_model(remove_orphans=True)            
        
        # Store the model properties in the array, write the model output, and increase the counter
        modelSizeDF.loc[curDir][0] = len(model.genes)
        modelSizeDF.loc[curDir][1] = len(model.metabolites)
        modelSizeDF.loc[curDir][2] = len(model.reactions)
        
    # Write the final model to a text file for inspection
    #    with open(modelDir+'/'+curDir+'/'+curDir+'Final.txt', 'w') as outFile:
    #        for curRxn in model.reactions:
    #            outFile.write(curRxn.id+'\t'+curRxn.build_reaction_string(use_metabolite_names=True)+'\n')
                
        print 'Processing model '+str(count)+' of '+str(len(dirList))
        cobra.io.write_sbml_model(model, modelDir+'/'+curDir+'/'+curDir+'.xml')
        count = count + 1
   
   # Write the results to file
    modelSizeDF.to_csv(summaryStatsDir+'/prunedModelStats.tsv', sep='\t')
     
    return
