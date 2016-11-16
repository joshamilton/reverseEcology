###############################################################################
# seedFunctions.py
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for working with seed sets.
################################################################################ 

# Import Python packages.
import os
import pandas as pd

# Define path for data included in the package
dataPath = os.path.dirname(os.path.abspath(__file__))+'/packageData'

################################################################################ 

# consolidateSeeds
# Consolidation seed weights into a single data frame.
# This function reads in the seed metabolites and their weights for each 
# genome. The data is read into a dataframe and then written to file. Structure:
# Rows: metabolites
# Columns: graphs
# Entries: unweighted seed set values

def consolidateSeeds(dirList, processedDataDir, summaryStatsDir):

    # Check that output directory exists
    if not os.path.exists(summaryStatsDir):
        os.makedirs(summaryStatsDir)

    print 'Consolidate seed sets'
        
    # Create a data frame to store the data
    revEcolMatrixDF = pd.DataFrame(columns=["Metabolite"])
    
    # Read in list of seed weights into a temporary data frame. Perform an outer
    # join with the existing data frame to incorporate the new list of weights.
    for curDir in dirList:
        tempDF = pd.read_csv(processedDataDir+'/'+curDir+'/'+curDir+'SeedCompounds.txt', names=['Metabolite','Name',curDir], sep='\t')
        tempDF =  tempDF.drop('Name', axis=1)
        revEcolMatrixDF = pd.merge(revEcolMatrixDF, tempDF, how='outer', on="Metabolite")
    
    # Replace all the NaN values with zeros
    revEcolMatrixDF.fillna(0, inplace=True)
    
    # Strip '_c' from metab names to facilitate the merge
    revEcolMatrixDF['Metabolite'] = revEcolMatrixDF['Metabolite'].str.rstrip('_c')
    
    # Append a new column containing common names associated with metabolite IDs.
    # The file metabMap.csv was created manually from the seed database, and should
    # be updated to reflect the particulars of your data set.
    namesDF = pd.read_csv(dataPath+'/metabMap.csv', names=['Metabolite','CommonName'])
    revEcolMatrixDF = pd.merge(revEcolMatrixDF, namesDF, how='inner', on="Metabolite")
    
    # Rearrange the order of the columns so that the common name is in front
    newOrder = revEcolMatrixDF.columns.tolist()
    newOrder = newOrder[-1:] + newOrder[:-1]
    revEcolMatrixDF = revEcolMatrixDF[newOrder]
    
    # Export the matrix of seed weights
    pd.DataFrame.to_csv(revEcolMatrixDF, summaryStatsDir+'/'+'seedMatrixWeighted.csv')    

    return revEcolMatrixDF
