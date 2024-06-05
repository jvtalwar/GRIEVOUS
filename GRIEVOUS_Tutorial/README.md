# GRIEVOUS Tutorial: Datasets and Validation!

Welcome intrepid investigator to the GRIEVOUS tutorial's datasets and validation directory. Here you will find the test datasets specified in the [tutorial](https://github.com/jvtalwar/GRIEVOUS/wiki/Data,-Preprocessing,-And-Formatting), in the `Datasets/Sarcoidosis/` subdirectory. You can download these datasets by cloning the GRIEVOUS github repository with `git clone https://github.com/jvtalwar/GRIEVOUS`. 

For validation of a complete tutorial `grievous` run, we provide two files:

    1) Tutorial_Intersecting_SNPs.tsv - A list of SNPs common to all 3 tutorial datasets, indexed according to CHROM:POS:sorted(REF:ALT) notation.
       - sorted(REF:ALT) corresponds to alphabetical order of nucleotides for a position. For example 1:33:A:C as opposed to 1:33:C:A.
    2) validate_tutorial.py - A script for comparing the equivalency of the output of your `grievous` tutorial run-through. Specifically this script takes in
                              the above Tutorial_Intersecting_SNPs.tsv file and the output of your `grievous intersect` run and assesses the equivalance in the outputs.
       - This script can also be used to assess the correctness of any given `grievous` run (post-intersect), provided a comprehensive set of intersecting SNPs in the 
         format of the Tutorial_Intersecting_SNPs.tsv file (this will need to be unique to your specific datasets of interest).


The run this validation script in accordance with the paths provided by the tutorial:

```bash
#PATH_TO_GRIEVOUS_TUTORIAL_INTERSECTING_SNPs=~/GRIEVOUS_Tutorial/Tutorial_Intersecting_SNPs.tsv #<-- Update this to the specific path of Tutorial_Intersecting_SNPs.tsv on your machine.
#PATH_TO_GRIEVOUS_INTERSECT_OUTPUT=~/GRIEVOUS_Tutorial/IntersectingVariants/AllIntersectingSNPs.tsv #<-- Update this with the specific path of your grievous tutorial intersect output. 
python validate_tutorial.py -i PATH_TO_GRIEVOUS_TUTORIAL_INTERSECTING_SNPs -g PATH_TO_GRIEVOUS_INTERSECT_OUTPUT
```
