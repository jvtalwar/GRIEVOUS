# @author: James V. Talwar

import pandas as pd
import logging
from collections import defaultdict
import os

logger = logging.getLogger(__name__)


'''
Converts a given file type to the GRIEVOUS Standard. Mapping file should at least contain all columns which deviate from standard.

Input(s):
   file: An unformatted GRIEVOUS column deviating pandas dataframe
   mappingFilePath: String corresponding to a mapping path 
Output(s):
    A formatted GRIEVOUS dataframe
'''
def FormatToStandard(file, mappingFilePath):
    mapping = pd.read_csv(mappingFilePath, delim_whitespace = True, header = None)
    mapping = dict(zip(mapping[0], mapping[1]))

    return file.rename(columns = mapping)

'''
About: This function creates a defacto graph where each CHR:POS can be considered a node and each ID an edge. 
       Non-singly connected nodes (i.e., nodes for which can reach another CHR:POS) are invalid and removed through 
       this preprocessing.

Input(s): Formatted pandas dataframe of a genomic object file (i.e., pvar or SSF)
Output(s): A processed dataframe of the genomic file filtered out for invalid SNPs that have multiple (RS)ID mappings
'''
def ExtractValidSNPs(formattedFile): 
    chrLocToID = defaultdict(set)
    idToChrLoc = defaultdict(set)
    
    for i in formattedFile.index: 
        chrLocToID[formattedFile.loc[i, "EasyCHRLOC"]].add(formattedFile.loc[i, "ID"]) 
        idToChrLoc[formattedFile.loc[i, "ID"]].add(formattedFile.loc[i, "EasyCHRLOC"]) 
        
    graph = defaultdict(lambda: defaultdict(set)) 
    
    for node,edges in chrLocToID.items():
        #Iterate through and populate the graph
        for edge in edges:
            graph[node][edge] = idToChrLoc[edge]
    
 
    keepForDownstreamProcessing = set()
    for node,edges in graph.items():
        validSNP = True
        for edge,otherNodes in edges.items():
            if not((len(otherNodes) == 1) and (next(iter(otherNodes)) == node)): #if there is ever a single instance of an ID pointing to more than one location then it is invalid and needs to be cut everywhere
                validSNP = False
        if validSNP:
            keepForDownstreamProcessing.add(node)
        
    toReturn = formattedFile[formattedFile.EasyCHRLOC.isin(keepForDownstreamProcessing)] #valid SNPs
    logger.info("{}/{} valid SNPs (with IDs not pointing to multiple locations) were extracted in preprocessing.\n".format(toReturn.shape[0], formattedFile.shape[0]))
        
    return toReturn

'''
Input(s): A dataframe of valid SNPs (i.e., SNPs where (RS)IDs are unique to a particular CHR LOC and are all valid SNPs in terms of REF/ALT) 
Output(s): A list of indexes of true biallelic SNPs 
'''
def ExtractBiallelicCandidates(chromPosDouble):
    logger.info("Beginning ExtractBiallelicCandidates() Filtering: {} SNPs under investigation... \n".format(chromPosDouble.shape[0]))
    
    #Step 0: Get sorted ref and alt allele --> This will be used in downstream filtering across file partitions
    chromPosDouble["SortedREFALT"] = ["".join(sorted(el)) for el in (chromPosDouble["REF"] + chromPosDouble["ALT"])]
    
    #Step 1: Are CHR, POS, REF, ALT jointly duplicated?
    amIDuplicatedByEverythingButID = chromPosDouble.duplicated(subset = ["CHR", "POS", "REF", "ALT"], keep = False)
    allDoubles = chromPosDouble[amIDuplicatedByEverythingButID]
    logger.info("Filter IDENTIFIED {} biallelic SNPs that were duplicated by CHR, POS, REF, ALT.".format(allDoubles.shape[0]))
    multiallelicSNPs = chromPosDouble[~amIDuplicatedByEverythingButID] #not really multiallelic SNPS, but not duplicated by chr,pos,ref,alt
    
    #Step 2: For those duplicated, everything matches but ID so drop duplicates:
    singlesTennis = allDoubles.drop_duplicates(subset = ["CHR", "POS", "REF", "ALT"]) #This can still have invalid SNPs such as double duplications
    
    #Step 3: Concatenate those SNPs that weren't jointly duplicated at all positions excepting ID with singlesTennis (i.e., those non-duplicated)
    stacked = pd.concat([singlesTennis, multiallelicSNPs], axis = 0)
    
    #Step 4: Partition into those duplicated by CHR and POS. Those unduplicated here are valid biallelic SNPs. 
    # Those duplicated need to continue to downstream filtering to identify any potential valid SNPs
    chrPosDuplication = stacked.duplicated(subset = ["CHR", "POS"], keep = False)
    noDups = stacked[~chrPosDuplication] #No duplicates by CHR, POS --> biallelic SNPs
    logger.info("Filter EXTRACTED {} biallelic SNPs that were unduplicated or were duplicated by CHR, POS, REF, ALT and were not duplicated through alternative means.".format(noDups.shape[0]))
    forFurtherFiltering = stacked[chrPosDuplication] #These are alleles that are duplicated by CHR, POS and need to proceed to step 5 
    
    del chrPosDuplication

    #Step 5: For those that were duplicated (dups) have to ask the question of are you unique by CHR, POS, SortedREFALT. If so then these are SNPs that were duplicated N number of times, but only with REF and ALT switched
    filteringDictionary = defaultdict(set) #Instantiate a checking dictionary
    for i in forFurtherFiltering.index:
        filteringDictionary[forFurtherFiltering.loc[i, "EasyCHRLOC"]].add(forFurtherFiltering.loc[i, "SortedREFALT"])
    
    realBiallelicChrLoc = set() 
    validIDs = {"AC", "AG", "AT", "CG", "CT", "GT"} #don't want any alleles that are of the form CC GG TT AA
    for k,v in filteringDictionary.items():
        if (len(v) == 1) and (next(iter(v)) in validIDs): #there is only one SortedREFALT orientation and its not of the doubled allele (e.g. CC) form
            realBiallelicChrLoc.add(k)
    
    #Get the remaining biallelic SNPs:
    sneakyBiallelics = forFurtherFiltering.EasyCHRLOC.isin(realBiallelicChrLoc)
    validSnpsAsWell = forFurtherFiltering[sneakyBiallelics] 
    logger.info("Filter IDENTIFIED AND EXTRACTED {} valid biallelic SNPs that were duplicated by CHR, POS due to duplication through REF/ALT switches.".format(validSnpsAsWell.shape[0]))

    invalidSnps = forFurtherFiltering[~sneakyBiallelics]  
    logger.info("Filter IDENTIFIED {} invalid SNPs that were duplicated by CHR, POS but had multiple alleles per position.".format(invalidSnps.shape[0]))
    
    #Now drop the duplicates from valid SNPs by those columns where duplication is relevant - namely CHR, POS, and sortedREFALT
    validSnpsAsWell = validSnpsAsWell.drop_duplicates(subset = ["CHR", "POS", "SortedREFALT"]) 
    
    
    #Step 6 (Nearly There, Nearly There...): Concatenate all the valid SNPs and return a list of indexes:
    toReturn = pd.concat([noDups, validSnpsAsWell], axis = 0)
     
    return list(toReturn.index)

'''
Input(s): A cleaned biallelic dataframe of a genomic file
Output(s): A dictionary mapping CHR:POS:REF:ALT to a given genomic file index
'''
def GenerateForwardIndices(genomicFile):
    goInThisOrder = ["CHR", "POS", "REF", "ALT"]
    grievousStandard = genomicFile[goInThisOrder].apply(lambda x: ":".join(map(str, x)), axis = 1)
    mapping = dict(zip(grievousStandard.values, grievousStandard.index))

    return mapping


'''
About: For all the keys that didn't map to grievous database/dictionary generate the alternative order - CHR:POS:ALT:REF

Input(s): temporalPincer: A dictionary of keys/indexes that were not found in grievous database (we're going both forward and backwards!)
Output(s): reversedGrievousKeys: A dictionary of CHR:POS:ALT:REF to the original input value (for which the key was CHR:POS:REF:ALT)
'''
def GenerateReverseIndices(temporalPincer):
    reversedGrievousKeys = defaultdict(str) 
    for k,v in temporalPincer.items():
        individualElements = k.split(":")
        if len(individualElements) != 4:
            logger.error("INDEX ERROR: Current index is length {} instead of 4. Exiting...".format(len(individualElements)))
            raise Exception("INDEX ERROR: Current index is length {} instead of 4. Exiting...".format(len(individualElements)))
            
        chrPosAltRef = ":".join([individualElements[0], individualElements[1], individualElements[3], individualElements[2]]) 
        reversedGrievousKeys[chrPosAltRef] = v #assign the reverse 
        
    return reversedGrievousKeys

'''
About: For all SNPs that are not indexed in a grievous chromosome db/dictionary in either orientation, ensure that the CHR:POS is not in said db/dictionary. If it is then the genomic SNP is invalid 
       as it points to different alleles than those of the db/dictionary. If not return a dictionary of CHR:LOC:REF:ALT: {CHR:STR POS:### REF:STR ALT:STR ID:STR} in the original genomic file 
       for that SNP.

Input(s):  1) snpsNotInDict: The biallelic genomic file SNPs not in the grievous db/dictionary (dataframe) 
           2) chromosomeDictionary: grievous chromosome-level db/dictionary (relating to the specific chromosome under investigation)
Output(s): 1) addendum: A dictionary of CHR:POS:REF:ALT --> {CHR:STR, POS:###, REF:STR, ALT:STR, ID:STR} for SNPs not in grievous db/dictionary 
           2) invalidSNPs: A set of SNP indexes where CHR POS is in the dictionary, but the alleles are different.
'''
def SanityCheckSNPsNotInDictionary(snpsNotInDict, chromosomeDictionary):
    goInThisOrder = ["CHR", "POS", "REF", "ALT"]
       
    uhOhSpudodios = snpsNotInDict[(snpsNotInDict.POS.isin(chromosomeDictionary.POS))] #Identified biallelic SNPs that are in the db/dictionary but with different SNPs at REF and ALT; chromosome-level so only need to filter by pos since all CHR match
    invalidSNPs = set(uhOhSpudodios.index) 
    validSNPsToAddToDict = snpsNotInDict[~snpsNotInDict.index.isin(uhOhSpudodios.index)]
    
    validSNPsToAddToDict.index = validSNPsToAddToDict[goInThisOrder].astype(str).apply(':'.join, axis=1)
    addendum = validSNPsToAddToDict[goInThisOrder + ["ID"]].T.to_dict() #subset down to only the 5 dictionary reqs  
    del validSNPsToAddToDict
    
    dictionaryCollisions = chromosomeDictionary[(chromosomeDictionary.POS.isin(uhOhSpudodios.POS))] 
    
    #Report in logs to user all SNPs that are identified as biallelic, but diverge from grievous database and thus are excluded
    for i in uhOhSpudodios.index:
        amIHere = dictionaryCollisions[dictionaryCollisions.POS == uhOhSpudodios.loc[i, "POS"]] 
        logger.info("SNP index {} (ID: {}) disagrees between pvar/ssf and database and will be removed. GRIEVOUS Database: {}; pvar/ssf: {}".format(i, uhOhSpudodios.loc[i, "ID"], list(amIHere.loc[:, "REF"] + amIHere.loc[:, "ALT"])[0], uhOhSpudodios.loc[i, "REF"] + uhOhSpudodios.loc[i, "ALT"]))
            
    del dictionaryCollisions
    
    return addendum, invalidSNPs

'''
About: Validates updated GRIEVOUS database/dictionary write contains all the biallelic SNPs not in the previous chromosome 
       alignment database/dictionary. 

Input(s):  1) updatedDictPath: String corresponding to the full path of the updated dictionary. 
           2) dictionaryAddendum: Dictionary corresponding to the biallelic SNPs (keys) and relevant information 
              that should have been added to the updated dictionary. 
           3) previousDict: Dataframe corresponding to the previous version (i.e., the one used for Orient()) of the 
              chromosome alignment dictionary.
           4) fileExtension: String in {parquet, tsv} corresponding to dictionary storage. Default parquet.
Output(s): validated: Boolean corresponding to the validation of the update dict 
'''
def ValidateDictUpdate(updatedDictPath, dictionaryAddendeum, previousDict, fileExtension = "parquet"): 
    assert fileExtension in {"parquet", "tsv"}, "Invalid database/dictionary fileExtension. Valid options are parquet or tsv."
    
    if fileExtension == "parquet":
        updatedDict = pd.read_parquet(updatedDictPath)
        updatedDict["CHR"] = updatedDict["CHR"].astype(str)
    else:
        updatedDict = pd.read_csv(updatedDictPath, sep = "\t", index_col = 0, dtype={'CHR': str})
            
    
    if len(updatedDict.index) != len(set(updatedDict.index)):
        logger.warning("Updated dictionary contains duplicated elements.")

    addOns = updatedDict[updatedDict.index.isin(dictionaryAddendeum.keys())]
    validated = addOns.shape[0] == len(dictionaryAddendeum)
    if not validated:
        logger.warning("Not all cohort identified biallelic SNPs were added to updated dictionary")
    
    validated = validated and (previousDict.shape[0] == (updatedDict.shape[0] - addOns.shape[0])) and (previousDict.shape[1] == updatedDict.shape[1])

    if not validated:
        logger.warning("Dimension inconsistency between previous and updated dictionaries.")

    return validated


'''
About: CleanReorientIndex is a method for handling duplications in pvar reorientation files. Specifically, it takes in a reorientIndex dataframe 
as extracted from Pvar.Write() and removes all duplicates (ensuring plink reorientation can be performed in cases of variant/index duplication).
The final output indexes can also be used to remove duplicates from the reorientRefAlleleThisWay dataframes. 
SSF file subset can also be used (as in SSF.Write()) to identify duplicates and generate Duplication_Reports.

Input(s):  1) duplicatedReorientIndex: reorientIndex dataframe (numerically indexed) unfiltered for duplicates.
           2) biallelicSNPs: list of identified chromosomal biallelic SNPs 
Output(s): 1) finalValidDF: A duplicate-free (i.e., filtered) reorientIndex dataframe
           2) duplicatedIDReport: dictionary corresponding to biallelic and non-biallelic IDs and their mapping to grievous index
           3) duplicatedIDCounts: dictionary corresponding to original file IDs and the number of times an ID occurs (if duplicated)
'''
def CleanReorientIndex(duplicatedReorientIndex, biallelicSnps):
    #Identify all duplicate IDs in the original file standard
    whereAreDuplications = duplicatedReorientIndex.duplicated(subset = ["ID"], keep = False)
    noDuplicates = duplicatedReorientIndex[~whereAreDuplications]
    checkForIssues = duplicatedReorientIndex[whereAreDuplications] 

    fastLookup = set(biallelicSnps) 

    #Generate a duplicated ID reports for the user:
    duplicatedIDReport = defaultdict(lambda: defaultdict(set))
    for i in checkForIssues.index:
        grievousIndex = checkForIssues.loc[i, "ProperlyFormattedName"]
        originalID = checkForIssues.loc[i, "ID"]
        if (grievousIndex in fastLookup):
            duplicatedIDReport["Biallelic_Duplicate"][originalID].add(grievousIndex)
        else:
            duplicatedIDReport["NOT_Biallelic_Duplicate"][originalID].add(grievousIndex)

    duplicatedCounts = checkForIssues.ID.value_counts()
    duplicatedIDCounts = dict(zip(duplicatedCounts.index, duplicatedCounts.values))
     
    noCompleteDuplicates = checkForIssues.drop_duplicates() #Drop duplicates if both original ID and grievous index are duplicated
    
    '''
     Now need to check if there are any duplicates by the original ID (i.e., first column) again - if so then there are biallelic SNPs doubled on the SNP chip and only 
     one is oriented the correct way as done during valid biallelic extraction to remove redundancies. For this to happen the original orientation would also need to 
     be opposite the dictionary orientation (if non-rsid; rsid duplication does not necessitate this).
    '''

    anyMoreDuplicates = noCompleteDuplicates.duplicated(subset = ["ID"], keep = False)
    canStackThese = noCompleteDuplicates[~anyMoreDuplicates]
    downstreamExtraction = noCompleteDuplicates[anyMoreDuplicates]
    
    canAlsoStackThese = pd.DataFrame()

    if downstreamExtraction.shape[0] != 0:
        logger.info(f"Duplication of original file ID, but unique grievous index case triggered. Identified {downstreamExtraction.shape[0]} instances\n")

        grievousBiallelicNumericalIndex = set()
        idHasBeenAdded = set()
    
        for i in downstreamExtraction.index: 
            grievousIndex = downstreamExtraction.loc[i, "ProperlyFormattedName"]
            if (grievousIndex in fastLookup):
                #elements are duplicated by ID, but for each duplicated ID each grievous index must be unique as have dropped complete duplicates
                grievousBiallelicNumericalIndex.add(i)

                # For instances of IDs duplicated one of which is marked as grievous biallelic (i.e., "ProperlyFormattedName" in biallelicSNPs) 
                # it is possible that the alternate orientation exists with the same ID. Need to identify those so they are not concatenated 
                # and thus duplications are included in NoDuplicates Files
                idHasBeenAdded.add(downstreamExtraction.loc[i, "ID"]) 
                

        biallelicCornerCase = downstreamExtraction[downstreamExtraction.index.isin(grievousBiallelicNumericalIndex)]

        # Handle remaining non-biallelic SNPs in downstream extraction without omission
        remainingDuplicates = downstreamExtraction[(~downstreamExtraction.index.isin(grievousBiallelicNumericalIndex)) & (~downstreamExtraction.ID.isin(idHasBeenAdded))].drop_duplicates(subset = ["ID"], keep = "first")

        canAlsoStackThese = pd.concat([biallelicCornerCase, remainingDuplicates], axis = 0)


    #stack, sort and return filtered reorient index
    finalValidDF = pd.concat([pd.concat([noDuplicates, canStackThese], axis = 0), canAlsoStackThese], axis = 0)
    finalValidDF = finalValidDF.sort_index()
    
    return finalValidDF, duplicatedIDReport, duplicatedIDCounts