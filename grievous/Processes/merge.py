#@author: James V. Talwar

import os
import pandas as pd
import logging
import re

logger = logging.getLogger(__name__)


'''
About: This method takes in a path to a grievous realigned write directory and generates composite reports across
all chromosomes that were completed at the time merge was called. If the write directory passed in was for an ssf
file, merge will also generate a composite ssf in the GRIEVOUS_Formatted subfolder.

Input(s):  1) realignment_directory: String corresponding to the --write_path (-w) of a completed grievous realign.
           2) updatedIDs: Boolean corresponding to whether to update/replace the original SSF IDs with the greivous realigned indexes.
Output(s): None
Write(s):  Composite Biallelic and Flipped SNP reports for a given cohort and a composite ssf if the realignmentDirectory was
           for an ssf file.
'''
def Merge(realignment_directory, updateIDs):
    #check if is a valid realignment directory
    assert os.path.isdir(realignment_directory), f"""Given path {realignment_directory} for --aligned_files does not exist/is not a directory. Please retry with a valid directory
    that is the output of a 'grievous realign' run."""

    subdirectories = os.listdir(realignment_directory)
    assert ("Reports" in subdirectories) and ('GRIEVOUS_Formatted' in subdirectories), f"Unable to locate both Reports and GRIEVOUS_Formatted folders in {realignment_directory}. Confirm directory - recall 'grievous realign' must be called before 'grievous merge'."

    #Check the file extension for ssfs - ensure unique (i.e., user did not reuse same directory for different realigns)
    grievousFormattedDir = os.path.join(realignment_directory, "GRIEVOUS_Formatted")
    fileType = set([el.split("_")[0].split("FormattedAndAligned")[1] for el in os.listdir(grievousFormattedDir) if el.endswith(".tsv")])

    assert (len(fileType) == 1) and ({"Pvar", "SSF"}.intersection(fileType) == fileType), f"{len(fileType)} file types found in --aligned_files directory {realignment_directory}. Given directory has realignment extensions: {fileType}. Ensure that every dataset-level grievous realign has a unique write_path."
    
    extension = next(iter(fileType))

    #grab all the biallelic SNPs and flipped SNPs and merge 
    reportDir = os.path.join(realignment_directory, "Reports")
    reportDirFiles = os.listdir(reportDir)

    if ("ALL_BiallelicSNPs.tsv" in reportDirFiles) or ("ALL_FlippedSNPs.tsv" in reportDirFiles):
        logger.warning(f"WARNING: Detected previous merge files ALL_BiallelicSNPs/ALL_FlippedSNPs.tsv in --aligned_files {os.path.join(realignment_directory, 'Reports')}. Overwriting...\n")
    
    mergeReportFiles = [os.path.join(reportDir, file) for file in reportDirFiles if (file.endswith(".tsv")) and (not re.match("WARNING_CHR.*_BiallelicDuplicates.tsv", file)) and (file not in {"ALL_FlippedSNPs.tsv", "ALL_BiallelicSNPs.tsv"})]

    biallelicSNPs = list()
    flippedSNPs = list()

    for report in mergeReportFiles:
        prefix = report.split("/")[-1].split("_")[-1].split(".")[0]
        assert prefix in {"BiallelicSNPs", "FlippedSNPs"}, f"Unknown report {report} with prefix {prefix} found. Valid prefixes for reports are 'BiallelicSNPs' and 'FlippedSNPs'. Exiting..."

        toAppend = pd.read_csv(report, sep = "\t")[prefix].tolist()

        if prefix == "BiallelicSNPs":
            biallelicSNPs = biallelicSNPs + toAppend

        else:
            flippedSNPs = flippedSNPs + toAppend


    #Write AllBiallelic and AllFlipped files
    if len(biallelicSNPs) > 0:
        pd.DataFrame(biallelicSNPs, columns = ["BiallelicSNPs"]).to_csv(os.path.join(reportDir, "ALL_BiallelicSNPs.tsv"), sep = "\t", index = False)

    if len(flippedSNPs) > 0:
        pd.DataFrame(flippedSNPs, columns = ["FlippedSNPs"]).to_csv(os.path.join(reportDir, "ALL_FlippedSNPs.tsv"), sep = "\t", index = False)

    logger.info("Biallelic and flipped SNP reports merged.")

    #if SSF merge chromosomal-level SSFs into one cohesive SSF
    if extension == "SSF":
        logger.info("\nMerging FormattedAndAlignedSSF files")
        mergeSSFFiles = [os.path.join(grievousFormattedDir, file) for file in os.listdir(grievousFormattedDir) if file.endswith(".tsv")]

        mergedSSF = pd.DataFrame()
        for ssf in mergeSSFFiles:
            chrSSF = pd.read_csv(ssf, sep = "\t", index_col = 0) 
            mergedSSF = pd.concat([mergedSSF, chrSSF], axis = 0)

        if updateIDs:
            mergedSSF.ID = mergedSSF.index
            mergedSSF.to_csv(os.path.join(grievousFormattedDir, "MergedSSF.ssf"), sep = "\t", index = False)
        
        else:
            mergedSSF.to_csv(os.path.join(grievousFormattedDir, "MergedSSF.ssf"), sep = "\t")

        
    logger.info("\ngrievous merge complete!")


    return None


