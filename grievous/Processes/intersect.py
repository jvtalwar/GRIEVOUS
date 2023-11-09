#@author: James V. Talwar

import os
import pandas as pd
import logging
import re
from .utils import _IdentifyPalindromic

logger = logging.getLogger(__name__)

'''
About: This method takes in n number of grievous realigned and merged cohorts and generates the feature/SNP set 
common to all cohorts of interest.

Input(s):  1) write_path: String corresponding to the directory to which write the grievous intersect feature set.
           2) cohorts: List of strings corresponding to the directories of grievous realigned and merged feature sets. 
           For each cohort this corresponds to the --write_path provided to grievous realign and to the --aligned_files 
           provided to grievous merge.
           3) exclude_biallelic_duplicates: Boolean corresponding to whether to exclude or include biallelic SNPs that were 
           originally duplicated in a cohort.
           4) omit_palindromic: Boolean corresponding to whether to omit palindromic SNPs from the intersected feature set.
           5) file_output_name: String corresponding to the output file name written by grievous intersect.    
Output(s): None
Write(s):  An intersected SNP/feature set common to all cohorts provided in directory write_path with file name
           file_output_name.
'''

def Intersect(write_path, cohorts, exclude_biallelic_duplicates, omit_palindromic, file_output_name):
    logger.info(f"Cohorts for intersection: {cohorts}\n")

    #Validate user inputs:
    assert os.path.isdir(write_path), f"--write_path {write_path} does not exist as a directory. Please retry with a valid write_path. Exiting..."
    
    for cohort in cohorts:
        assert os.path.isdir(cohort), f"Provided cohort {cohort} does not exist as a directory. Exiting..."

        subdirectories = os.listdir(cohort)
        assert ("Reports" in subdirectories) and ('GRIEVOUS_Formatted' in subdirectories), f"""Unable to locate both Reports and GRIEVOUS_Formatted folders in cohort {cohort}.
          Confirm directory and recall 'grievous intersect' is called after a given cohort has undergone both 'grievous realign' and 'grievous merge'."""
        
        assert "ALL_BiallelicSNPs.tsv" in os.listdir(os.path.join(cohort, "Reports")), f"""ALL_BiallelicSNPs.tsv file missing from cohort {cohort}. Validate 'grievous merge' was
        run. If so, then no biallelic SNPs exist for cohort {cohort} and the intersection across all cohorts of interest is an empty set."""

    intersectingSNPs = set()
    firstInstance = True

    for cohort in cohorts:
        reportDir = os.path.join(cohort, "Reports")
        cohortAlignedBiallelicSNPs = set(pd.read_csv(os.path.join(reportDir, "ALL_BiallelicSNPs.tsv"), sep = "\t")["BiallelicSNPs"])

        if exclude_biallelic_duplicates:
            duplicatesToOmit = set()
            cohortBiallelicDuplicateFiles = [os.path.join(reportDir, file) for file in os.listdir(reportDir) if re.match("WARNING_CHR.*_BiallelicDuplicates.tsv", file)]
            for dups in cohortBiallelicDuplicateFiles:
                 duplicatesToOmit = duplicatesToOmit.union(set(pd.read_csv(dups, sep = "\t", header = None)[0]))

            cohortAlignedBiallelicSNPs = cohortAlignedBiallelicSNPs.difference(duplicatesToOmit)

        if firstInstance: # base case (i.e., first cohort in cohorts)
            firstInstance = False
            
            #Handle corner case where all SNPs are duplicated and exclude_biallelic_duplicates called - write an empty intersection (no SNPs to intersect)
            if len(cohortAlignedBiallelicSNPs) == 0:
                break

            intersectingSNPs = cohortAlignedBiallelicSNPs

        else:
            intersectingSNPs = intersectingSNPs.intersection(cohortAlignedBiallelicSNPs)


    #Handle palindromic SNP exclusion:
    if omit_palindromic:
        logger.info("\nExcluding palindromic SNPs...\n")
        palindromicSNPs = _IdentifyPalindromic(intersectingSNPs)
        intersectingSNPs = intersectingSNPs.difference(palindromicSNPs)
        
        logger.info(f"\nIdentified and removed {len(palindromicSNPs)} palindromic SNPs.\n")

    logger.info(f"The number of intersecting SNPs across the given {len(cohorts)} cohorts is: {len(intersectingSNPs)}. Writing intersecting SNPs...")

    pd.DataFrame(intersectingSNPs).to_csv(os.path.join(write_path, file_output_name + ".tsv"), sep = "\t", index = False, header = None)

    logger.info(f"\nIntersecting SNPs/features identified and written to {os.path.join(write_path, file_output_name + '.tsv')}. My work here is done...")

    return None