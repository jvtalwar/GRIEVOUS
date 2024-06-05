import pandas as pd
import logging
import argparse

'''
About: Method to validate grievous' correctness on tutorial data.

Input(s): 1) tutorial_intersecting_snps: String corresponding to the path to all intersecting features (identified independently of grievous).
                                         These features should be stored as a .tsv file, with IDs corresponding to CHR:POS:sorted(REF:ALT) notation.
                                         sorted(REF:ALT) corresponds to alphabetical order of nucleotides for a position. For example 1:33:A:C as opposed
                                         to 1:33:C:A.
          2) grievous_intersecting_snps: String corresponding to the file path of output of grievous intersect for the tutorial (or equivalent) data. 
                                         If validating grievous with independent dataset ensure tutorial_intersecting_snps follows the aforementioned 
                                         standards. 
                                    
'''
def Validate_Tutorial(tutorial_intersecting_snps, grievous_intersecting_snps):
    generatedIntersectingFeatures = set(pd.read_csv(tutorial_intersecting_snps, sep = "\t", header = None)[0].tolist())
    grievousIntersect = set(pd.read_csv(grievous_intersecting_snps, sep = "\t", header = None)[0].tolist())

    #Ensure tutorialIntersectingFeatures are stored in nucleotide alphabetical order:
    for snp in generatedIntersectingFeatures:
        chrom, pos, ref, alt = snp.split(":")
        assert [ref, alt] == sorted(ref + alt), "tutorial intersecting snps are not stored in nucleotide alphabetical order."

    
    #convert grievousIntersect snps to alphabetical order. Recall grievous orients new SNPs to a dataset orientation (not alphabetically).
    sortedGrievousIntersect = set()
    for snp in grievousIntersect:
        chrom, pos, ref, alt = snp.split(":")
        sortedRef, sortedAlt = sorted(ref + alt)
        sortedGrievousIntersect.add(":".join([chrom, pos, sortedRef, sortedAlt]))

    return sortedGrievousIntersect == generatedIntersectingFeatures, len(sortedGrievousIntersect) == len(generatedIntersectingFeatures), len(generatedIntersectingFeatures), len(sortedGrievousIntersect)


def main():
    logger = logging.getLogger()
    handler = logging.StreamHandler()
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    parser = argparse.ArgumentParser()

    path_args = parser.add_argument_group("Input options:")
    path_args.add_argument('--independent_intersecting_snps', '-i', type = str, help= 'Path to (tutorial) cross-dataset intersecting snps identified independently of grievous.')
    path_args.add_argument("--grievous_intersecting_snps", '-g', type = str, help = "Path to the output of grievous intersect for all datasets employed to identify independent_intersecting_snps.")
    
    args = parser.parse_args()

    equivalent, sameNumberIdentified, num_snps_ind, num_snps_grievous = Validate_Tutorial(tutorial_intersecting_snps = args.independent_intersecting_snps,
                                   grievous_intersecting_snps = args.grievous_intersecting_snps)
    

    if equivalent:
        logger.info("\nPASS!!!\ngrievous intersect output matches independently identified intersecting snps.\n")

    else:
        logger.warning("\nFAIL!!!\nindependent_intersecting_snps and grievous_intersecting_snps are not equivalent.\n")

        if not sameNumberIdentified:
            logger.warning(f"Number of GRIEVOUS intersect identified biallelic snps across datasets: {num_snps_grievous}")
            logger.warning(f"Number of provided independent intersecting biallelic snps across datasets: {num_snps_ind}")

if __name__ == "__main__":
    main()