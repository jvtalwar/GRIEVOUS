# @author: James V. Talwar

import pandas as pd
import os
import logging
from ._utils import (FormatToStandard, ExtractValidSNPs, ExtractBiallelicCandidates, GenerateForwardIndices, 
                     GenerateReverseIndices, SanityCheckSNPsNotInDictionary, ValidateDictUpdate, CleanReorientIndex)

logger = logging.getLogger(__name__)

class Pvar:
    '''
    Parameters:
       filePath: String corresponding to where the file can be found (e.g., pvar)
       mappingPath: String corresponding to where column conversion standard can be found (default = None)
       commentCharsToSkip: List of strings corresponding to contiguous leading comments as are found across genomic files (default = ["##"]). 
                           
    '''
    def __init__(self, filePath, mappingPath = None, commentCharsToSkip = ["##"], *args, **kwargs):  
        #Handle when the top of the file is commented
        numLinesToSkip = 0
        with open(filePath, "r") as f:
            for line in f:
                isNotCommented = True
                for commentChar in commentCharsToSkip:
                    if line[:len(commentChar)] == commentChar:
                        numLinesToSkip += 1
                        isNotCommented = False
                        break
                
                if isNotCommented:
                    break

        if mappingPath:
            file = pd.read_csv(filePath, sep = "\t", skiprows = numLinesToSkip)
            self.file = FormatToStandard(file, mappingPath)
        else:
            self.file = pd.read_csv(filePath, sep = "\t", skiprows = numLinesToSkip)
        
        #Ensure the all key columns exist and are properly formatted to GRIEVOUS standard
        assert {"CHR", "POS", "ID", "REF", "ALT"}.intersection(set(self.file.columns)) == {"CHR", "POS", "ID", "REF", "ALT"}, "File formatting error: Mandatory columns [CHR, POS, ID, REF, ALT] are missing. Exiting..."

        self.file["CHR"] = self.file["CHR"].astype(str) #ensure consistency when handling X chromosome
        assert len(set(self.file.CHR)) == 1, "More than 1 Chromosome found in given file. Alignment expects files at the chromosome resolution. Exiting..."
        assert self.file.CHR[0] in ([str(i) for i in range(1,23)] + ["X", "Y", "MT"]), f"File CHR {self.file.CHR[0]} deviates from grievous standard. Exiting..."

        self.fileType = ".pvar" #define file type for inheritance/reusability purposes

        #Define internal class states which will be updated during the GRIEVOUS alignment process
        self._candidateBiallelicSNPIDs = None #Candidate biallelic index subset of the genomic file as updated by the Clean() method below
        self._theseIndexesNeedToFlipTheirRefAndAltAlleles = None #Indexes of SNPs that need to be reoriented - updated by Orient() below
        self._biallelicSNPIndexes = None #The GRIEVOUS identified biallelic index SNP indexes (post QC and dictionary divergence check)
        
        #Key states:  
        self.dictionaryAddendum = dict() #updated in Orient()

        #update in Align()
        self.biallelicSNPs = list()
        self.flippedSNPs = list()

    @property
    def Chrom(self):
        return self.file.CHR[0]
    
    
    '''
    About: Perform data quality filtering and generate a biallelic candidate SNP subset. 

    Input(s): None 
    Output(s): None
    Update(s): 1) self._candidateBiallelicSNPIDs: List of self.file indexes that pass genomic file only (i.e., not dictionary compared) biallelic checks.         

    '''
    def Clean(self):
        logger.info("Total number of SNPs in this Chromosome {}: {}\n".format(self.fileType, self.file.shape[0]))

        self.file["EasyCHRLOC"] = [el for el in (self.file["CHR"] + ":" + self.file["POS"].astype(str))] # helper column for ease of downstream Clean method processing
        
        #Ensure valid SNPs --> Remove all those that are not in biallelic vocabulary of {A,C,G,T}; This will also remove all variants that are not length of 1
        validSNPs = {"A", "C", "G", "T"}
        validFile = self.file[(self.file.REF.isin(validSNPs)) & (self.file.ALT.isin(validSNPs))]

        #Preprocess - remove all SNPs with IDs pointing to multiple locations
        idCleanedFile = ExtractValidSNPs(validFile)

        #Perform QC on idCleanedFile to get a set of biallelic SNP candidates (which will be compared against the grievous db/dictionary in Orient)
        biallelicCandidateIndices = ExtractBiallelicCandidates(idCleanedFile)
        
        #Check for potential AA, CC, TT, GG, alleles and remove
        ensureThingsAreValid = validFile.loc[biallelicCandidateIndices, :]
        sanityChecker = ensureThingsAreValid[ensureThingsAreValid.REF == ensureThingsAreValid.ALT].shape[0] == 0
        logger.info("Sanity Check: Are all identified biallelic candidate SNPs valid (i.e., not AA, CC, GG, TT)? {}\n".format(sanityChecker))
        
        if not sanityChecker:
            logger.info("{} AA, CC, GG, TT SNPs detected. Removing now...\n".format(ensureThingsAreValid[ensureThingsAreValid.REF == ensureThingsAreValid.ALT].shape[0]))
            ensureThingsAreValid = ensureThingsAreValid[ensureThingsAreValid.REF != ensureThingsAreValid.ALT]
    
        self._candidateBiallelicSNPIDs = list(ensureThingsAreValid.index)
        self.file = self.file.drop(["EasyCHRLOC"], axis = 1) #remove helper column

        return None

    '''
    About: Orient SNPs to the GRIEVOUS database/dictionary standard. Identify those SNPs which exist in both the db/dictionary orientation (and the reverse orientation). 
    Perform final QC on biallelic SNPs that are not in db/dictionary before adding them as biallelic SNPs and to future db/dictionary versions.
    
    Input(s): 1) chromosomeDictionary: chromosome level dataframe 
    Output(s): None
            
    Update(s):
            1) self._theseIndexesNeedToFlipTheirRefAndAltAlleles: List of given cohort genomic file indexes which exist in the opposite orientation from db/dictionary 
            2) self._biallelicSNPIndexes: List of all biallelic SNP indexes as identified by GRIEVOUS
            3) self.dictionaryAddendum: Dictionary of SNPs to be added to future GRIEVOUS db/dictionary versions

    '''
    def Orient(self, chromosomeDictionary):
        biallelicCandidateSubset = self.file.loc[self._candidateBiallelicSNPIDs, :] #true bialleleic SNPs
        
        #convert the dictionary indexes to a set for faster lookup
        doIExist = set(chromosomeDictionary.index)
        grievousOrientationToIndex = GenerateForwardIndices(biallelicCandidateSubset) 

        #Ask the question of which elements are not found in the grievous db/dictionary
        problemChildren = set([k for k in grievousOrientationToIndex]).difference(doIExist) #all the CHR:LOC:REF:ALT positions that don't match the grievous db/dictionary 
        logger.info("{}/{} SNPs were found in grievous chromosome database in their given orientation.".format(len(grievousOrientationToIndex) - len(problemChildren), biallelicCandidateSubset.shape[0]))

        #Generate keys in the other orientation
        twistIt = GenerateReverseIndices({k:grievousOrientationToIndex[k] for k in problemChildren})
        
        #Find which keys are in the dictionary but in the other orientation (and thus need to be flipped in the pvar file)
        reverseOrientationInDict = set([k for k in twistIt]).intersection(doIExist)
        self._theseIndexesNeedToFlipTheirRefAndAltAlleles = [twistIt[k] for k in reverseOrientationInDict] #These are the indexes of SNPs that need to be changed 
        logger.info("{}/{} SNPs were found in grievous chromosome database in the opposite orientation.".format(len(self._theseIndexesNeedToFlipTheirRefAndAltAlleles), biallelicCandidateSubset.shape[0]))
        
        #Get the remaining SNPs (i.e., those that are not found in mega dictionary in either orientation)
        newChallengersApproach = set([k for k in twistIt]).difference(doIExist)
        indexOfSnpsThatNeedToBeSanityCheckedBeforeAddingToMegaDictionary = [twistIt[k] for k in newChallengersApproach]
        logger.info("{}/{} SNPs were not found in grievous chromosome database. These SNPs are proceeding to db/dictionary addition validation.\n".format(len(indexOfSnpsThatNeedToBeSanityCheckedBeforeAddingToMegaDictionary), biallelicCandidateSubset.shape[0]))

        dictionaryAddendum, snpIndexesThatDisagreeWithDictionary = SanityCheckSNPsNotInDictionary(biallelicCandidateSubset.loc[indexOfSnpsThatNeedToBeSanityCheckedBeforeAddingToMegaDictionary, :], chromosomeDictionary) 
        logger.info("There were {} genomic biallelic SNPs that deviated from grievous chromosome database. Removing these...\n".format(len(snpIndexesThatDisagreeWithDictionary)))

        #Get the list of valid SNP indexes/remove invalid SNPs (those that were in mega dict, but differed by allele)
        self._biallelicSNPIndexes = [i for i in biallelicCandidateSubset.index if i not in snpIndexesThatDisagreeWithDictionary]
        
        #Update internal dictionary state
        self.dictionaryAddendum = dictionaryAddendum
          
        return None

    '''
    About: GRIEVOUS database/dictionary ALIGNs (true biallelic SNPs incorrectly oriented are swapped) the full genomic file (self.file) indexed by CHR:POS:REF:ALT (after correction). 
    It also extracts and updates the true biallelic indexes by CHR:LOC:REF:ALT (converting from the numerical indexes). 
     - NOTA BENE: This method should be called after Orient() which updates the necessary internal states self._theseIndexesNeedToFlipTheirRefAndAltAlleles and self._biallelicSNPIndexes.

    Input(s):  1) writePath: String corresponding to the directory for which to write duplication warning files 
    Output(s): None

    Update(s): 1) self.file: Aligned/Corrected and properly indexed file
               2) self.biallelicSNPs: A list of CHR:LOC:REF:ALT biallelic SNP indexes used for merging to full cohort biallelic SNP

    Write(s):  1) Warning file if any true biallelic SNPs exist as duplicated in the underlying genomic file 
'''
    def Align(self, writePath):
        assert self._biallelicSNPIndexes is not None, ".Orient() must be called before .Align()"

        #Swap alleles that are backwards and orient the Beta coefficient correctly if differs from db/dictionary
        chromosomeNumber = self.file.CHR[0] 

        swapCondition = self.file.index.isin(self._theseIndexesNeedToFlipTheirRefAndAltAlleles)
        self.file.loc[swapCondition, ["REF", "ALT"]] = self.file.loc[swapCondition, ["ALT", "REF"]].values

        if self.fileType == ".ssf": #SSF files (from which this class inherits) need to also reorient their Betas when in the opposite orientation
            self.file.loc[swapCondition, "BETA"] *= -1

        #Index the dataframe according to CHR:POS:REF:ALT now that it has been db/dictionary oriented
        goInThisOrder = ["CHR", "POS", "REF", "ALT"]
        self.file.index = [":".join([str(self.file.loc[i,col]) for col in goInThisOrder]) for i in self.file.index]
        
        #Now that indexed, oriented, and aligned get the biallelic CHR:POS:REF:ALT indexes
        formattedBiallelicIndexes = list(map(list(self.file.index).__getitem__, self._biallelicSNPIndexes))
        flippedSNPs = list(map(list(self.file.index).__getitem__, self._theseIndexesNeedToFlipTheirRefAndAltAlleles))

        #Sanity check: This can happen when SNPs are completely duplicated in the original pvar by CHR POS REF ALT; Log/Report warning SNPs to user
        if self.file[self.file.index.duplicated(keep = False)].shape[0] != 0: 
            whoseCausingProblems = self.file[self.file.index.duplicated(keep = False)]
            logger.warning("WARNING: {} SNP duplication events exist in the original file.".format(whoseCausingProblems.shape[0])) 
            biallelicProblems = set(whoseCausingProblems.index).intersection(set(formattedBiallelicIndexes))
            logger.warning("Of the SNPs duplicated by GRIEVOUS formatted index, {} of them is/are biallelic SNPs.\n".format(len(biallelicProblems)))
            if len(biallelicProblems) > 0:
                logger.warning("These biallelic SNPs are {}. Writing to WARNING_CHR{}_BiallelicDuplicates file now.\n".format(biallelicProblems, chromosomeNumber))
                pd.DataFrame(biallelicProblems).to_csv(os.path.join(writePath, "Reports/WARNING_CHR{}_BiallelicDuplicates.tsv".format(chromosomeNumber)), sep = "\t", index = False, header = None)
            #logger.warning("All duplicated index elements are:\n{}".format(whoseCausingProblems.index))
            

        #Update object biallelic SNPs and flipped SNPs
        self.biallelicSNPs = formattedBiallelicIndexes
        self.flippedSNPs = flippedSNPs

        return None 
        
    '''
    About: Updates a GRIEVOUS chromosome alignment database/dictionary after cohort orientation by writing the updated
    database/dictionary, validating the update, and then replacing the previous chromosome database/dictionary with the updated one.  

    Input(s): 1) chromosomeDictionary: Dataframe corresponding to the chromosome dictionary used for Orient().
              2) dictPath: String corresponding to the directory where the chromosome dictionary is kept.
              3) fileExtension: String in {parquet, tsv} corresponding to dictionary storage. Default parquet.
    Output(s): None
    Write(s): 
              1) An updated chromosome alignment dictionary with all biallelic SNPs for the given chromosome (identified in a 
                 given cohort, and thus ensuring orientation consistency for all subsequent cohorts that have these SNPs as well).
    '''
    def UpdateDict(self, chromosomeDictionary, dictPath, fileExtension = "parquet"):
        assert self._biallelicSNPIndexes is not None, ".Orient() must be called before .UpdateDict()"
        assert fileExtension in {"parquet", "tsv"}, "Invalid database/dictionary fileExtension. Valid options are parquet or tsv."
        assert f"CHR_{self.Chrom}.{fileExtension}" in os.listdir(dictPath), "Chromosome alignment database/dictionary not found at dictPath."

        if len(self.dictionaryAddendum) == 0:
            logger.info(f"All biallelic SNPs in CHR {self.Chrom} are in the current alignment database/dictionary. Update is not required.")
            return None

        updateName = f"CHR_{self.Chrom}_UPDATED.{fileExtension}"
        updatedDictPath = os.path.join(dictPath, updateName)

        #update the dictionary and write:
        logger.info("Adding {} total SNPs to GRIEVOUS database".format(len(self.dictionaryAddendum)))
        updatedDict = pd.concat([chromosomeDictionary, pd.DataFrame(self.dictionaryAddendum).T], axis = 0)

        if fileExtension == "parquet":
            updatedDict.to_parquet(updatedDictPath)
        else:
            updatedDict.to_csv(updatedDictPath, sep = "\t")

        #Validate the updated dictionary and replace the previous version with the updated version
        validated = ValidateDictUpdate(updatedDictPath, self.dictionaryAddendum, chromosomeDictionary, fileExtension)

        if not validated:
            os.remove(updatedDictPath)
            raise Exception(f"ERROR: All biallelic SNPs not found in previous CHR {self.Chrom} alignment dictionary were not written successfully. Exiting...")
        
        else:
            logger.info("Dictionary update validated! Replacing previous alignment dictionary with the updated version.\n")
            os.remove(os.path.join(dictPath, f"CHR_{self.Chrom}.{fileExtension}"))
            os.rename(updatedDictPath, os.path.join(dictPath, f"CHR_{self.Chrom}.{fileExtension}"))

        return None
    

    '''
    About: Write all reports as well as GRIEVOUS Formatted and reorientation files, the latter of which can be passed to 
           PLINK2 to align the underlying plink genotype files. Reorientation files are exclusive to self.fileType == '.pvar'  
    
    Input(s): 1) writePath: String corresponding to the directory where want to write these files.
              2) returnFullFile: Boolean corresponding to whether to return the full aligned pvar file (including
              any uncorrected information such as INFO columns) or the grievous PVAR column subset: [CHR, POS, ID, REF, ALT].
    
    Output(s): None

    Write(s): 
        1) The chromosome-level aligned pvar 
        2) The chromosome-level ref-allele files  (can be passed to PLINK2 after duplicate removal)
        3) The chromosome-level original pvar ID to formatted pvar index files (can be passed to PLINK2 after duplicate removal)
        4) A chromosome-level original ID SNP duplication report (if any duplicate IDs exist in the original file)
        5) The chromosome-level identified biallelic SNPs
        6) The chromosome-level subset of identified biallelic SNPs that were reoriented from cohort orientation to grievous dictionary standard (if these exist).
        
    '''
    def Write(self, writePath, returnFullFile):
        getChromosomeID = self.file.CHR[0] 
        
        #Write the formatted chromosome-level aligned pvar file (can be used to ensure consistency with new files generated by plink)
        if returnFullFile:
            self.file.to_csv(os.path.join(writePath, "GRIEVOUS_Formatted/FormattedAndAlignedPvar_CHR{}.tsv".format(getChromosomeID)), sep = "\t")
        else:
            self.file[["CHR", "POS", "ID", "REF", "ALT"]].to_csv(os.path.join(writePath, "GRIEVOUS_Formatted/FormattedAndAlignedPvar_CHR{}.tsv".format(getChromosomeID)), sep = "\t")

        #create and write reorientation files (including duplicates)
        logger.info("Writing PLINK2 update files... In total there are {} SNPs for chromosome {}.".format(self.file.shape[0], getChromosomeID))

        refAlleleFile = self.file[["ID", "REF"]]
        refAlleleFile.to_csv(os.path.join(writePath, "Reorientation/ReorientRefAlleleThisWay_CHR{}.tsv".format(getChromosomeID)), sep = "\t", index = False, header = None)

        logger.info("Finished writing ref allele file...")
    
        oldIDToFormattedID = self.file[["ID"]].copy()
        oldIDToFormattedID["ProperlyFormattedName"] = oldIDToFormattedID.index
        oldIDToFormattedID.to_csv(os.path.join(writePath, "Reorientation/ReorientIndex_CHR{}.tsv".format(getChromosomeID)), sep = "\t", index = False, header = None)

        logger.info("Finished writing ID update file...")

        #Generate cleaned non-duplicated ReorientRefAlleleThisWay_ and ReorientIndex_ files for Plink realignment
        oldIDToFormattedID.index = [i for i in range(len(oldIDToFormattedID.index))]
        noDuplicateReorientIndex, duplicatedIDReport, duplicatedIDCounts = CleanReorientIndex(oldIDToFormattedID, self.biallelicSNPs)

        refAlleleFile.index = [i for i in range(len(refAlleleFile.index))]
        noDuplicateRefAlleleFile = refAlleleFile[refAlleleFile.index.isin(noDuplicateReorientIndex.index)] #numerical indexing is the same across the 2 files as both subset from self.file columns

        noDuplicateReorientIndex.to_csv(os.path.join(writePath, "Reorientation/NoDuplicates_ReorientIndex_CHR{}.tsv".format(getChromosomeID)), sep = "\t", index = False, header = None)
        noDuplicateRefAlleleFile.to_csv(os.path.join(writePath, "Reorientation/NoDuplicates_ReorientRefAlleleThisWay_CHR{}.tsv".format(getChromosomeID)), sep = "\t", index = False, header = None)
        
        #Generate duplication report (if there are indeed duplicates):
        if len(duplicatedIDCounts) > 0:
            with open(os.path.join(writePath, "Reports/CHR{}_Duplication_Report.txt".format(getChromosomeID)), "a") as duplicationReport:
                duplicationReport.write("Original file duplication ID incidences:\n\n")
                
                for k,v in duplicatedIDCounts.items():
                    duplicationReport.write("\t".join([k,str(v)]) + "\n")
                
                duplicationReport.write("\nOriginal ID to grievous index unique incidence mappings:\n")

                for k,v in duplicatedIDReport.items():
                    duplicationReport.write("\n" + k + "\n\n")
                    for originalID, grievousIndexes in v.items():
                        duplicationReport.write("\t".join([originalID , str(grievousIndexes)]) + "\n")   


        logger.info("Duplicate removal for ReorientIndex and ReorientRefAlleleThisWay complete...\n")

        logger.info("Writing biallelic SNPs and subset of biallelic SNPs that were reoriented from cohort specific orientation to grievous database standard\n")
        
        pd.DataFrame(self.biallelicSNPs, columns = ["BiallelicSNPs"]).to_csv(os.path.join(writePath, "Reports/CHR{}_BiallelicSNPs.tsv".format(getChromosomeID)), sep = "\t", index = False)
        if len(self.flippedSNPs) > 0:
            pd.DataFrame(self.flippedSNPs, columns = ["FlippedSNPs"]).to_csv(os.path.join(writePath, "Reports/CHR{}_FlippedSNPs.tsv".format(getChromosomeID)), sep = "\t", index = False)

        logger.info("GRIEVOUS Alignment for CHR {} complete\n".format(getChromosomeID))

        return None