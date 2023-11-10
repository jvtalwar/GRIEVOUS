# @author: James V. Talwar

import pandas as pd
import os
import logging
from .pvar import Pvar
from ._utils import CleanReorientIndex

logger = logging.getLogger(__name__)

class SSF(Pvar):
    '''
    Parameters:
       filePath: String corresponding to where the file can be found (e.g., ssf)
       mappingPath: String corresponding to where column conversion standard can be found (default = None)
       commentCharsToSkip: List of strings corresponding to contiguous leading comments as are found across genomic files (default = ["##"]). 
                           
    '''
    def __init__(self, filePath, mappingPath = None, commentCharsToSkip = ["##"], *args, **kwargs):
        super().__init__(filePath, mappingPath, commentCharsToSkip, *args, **kwargs) 
        self.fileType = ".ssf"

        assert "BETA" in self.file.columns, "SSF Beta value missing (expecting BETA)"

    '''
    About: Writes all reports as well as generates a CHR-level formatted, oriented, and aligned SSF file needed for tasks such as SNP-set identification 
           and extraction by p-value in genotype cohorts or PRS construction. Note this CHR-level SSF employs the aligned GRIEVOUS standard IDs as the
           index and returns the original IDs as ID.
    
    Input(s):  1) writePath: String corresponding to the directory where want to write file(s).
               2) returnFullFile: Boolean corresponding to whether to return the full aligned ssf file (including
                  any uncorrected information such as INFO columns) or the grievous SSF column subset: [CHR, POS, ID, REF, ALT, BETA].
    
    Output(s): None 
    
    Write(s): 1) Chromosome-level grievous formatted and aligned SSF (as a .tsv) 
              2) A chromosome-level original ID SNP duplication report (if any duplicate IDs exist in the original file)
              3) The chromosome-level identified biallelic SNPs
              4) The chromosome-level subset of identified biallelic SNPs that were reoriented from cohort orientation to grievous dictionary standard (if these exist).
    ''' 
    def Write(self, writePath, returnFullFile):
        getChromosomeID = list(self.file.CHR)[0]
        
        if returnFullFile:
            self.file.to_csv(os.path.join(writePath, "GRIEVOUS_Formatted/FormattedAndAlignedSSF_CHR{}.tsv".format(getChromosomeID)), sep = "\t")
        
        else:
            self.file[["CHR", "POS", "ID", "REF", "ALT", "BETA"]].to_csv(os.path.join(writePath, "GRIEVOUS_Formatted/FormattedAndAlignedSSF_CHR{}.tsv".format(getChromosomeID)), sep = "\t")
    
        logger.info("Finished writing aligned and formatted chromosome {} ssf\n".format(getChromosomeID))

        #Identify and report ID duplicates in the SSF if they exist
        if sum(self.file.duplicated(subset = ["ID"], keep = False)) > 0:
            logger.info(f"Original CHR {self.Chrom} SSF file has duplicate IDs - generating duplication report.")

            oldIDToFormattedID = self.file[["ID"]].copy()
            oldIDToFormattedID["ProperlyFormattedName"] = oldIDToFormattedID.index
            oldIDToFormattedID.index = [i for i in range(len(oldIDToFormattedID.index))]

            _, duplicatedIDReport, duplicatedIDCounts = CleanReorientIndex(oldIDToFormattedID, self.biallelicSNPs)
            
            with open(os.path.join(writePath, "Reports/CHR{}_Duplication_Report.txt".format(getChromosomeID)), "a") as duplicationReport:
                duplicationReport.write("Original file duplication ID incidences:\n\n")
                
                for k,v in duplicatedIDCounts.items():
                    duplicationReport.write("\t".join([k,str(v)]) + "\n")
                
                duplicationReport.write("\nOriginal ID to grievous index unique incidence mappings:\n")

                for k,v in duplicatedIDReport.items():
                    duplicationReport.write("\n" + k + "\n\n")
                    for originalID, grievousIndexes in v.items():
                        duplicationReport.write("\t".join([originalID , str(grievousIndexes)]) + "\n")

        logger.info("Writing biallelic SNPs and subset of biallelic SNPs that were reoriented from cohort specific orientation to grievous database standard\n")
        
        pd.DataFrame(self.biallelicSNPs, columns = ["BiallelicSNPs"]).to_csv(os.path.join(writePath, "Reports/CHR{}_BiallelicSNPs.tsv".format(getChromosomeID)), sep = "\t", index = False)
        if len(self.flippedSNPs) > 0:
            pd.DataFrame(self.flippedSNPs, columns = ["FlippedSNPs"]).to_csv(os.path.join(writePath, "Reports/CHR{}_FlippedSNPs.tsv".format(getChromosomeID)), sep = "\t", index = False)

        logger.info("GRIEVOUS Alignment for CHR {} complete\n".format(getChromosomeID))
        
        return None
