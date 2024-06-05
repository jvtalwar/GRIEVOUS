# GRIEVOUS Tutorial: Datasets and Validation

Welcome intrepid investigator to the GRIEVOUS tutorial's datasets and validation directory. 

## Downloading GRIEVOUS Tutorial Datasets:

To download the (toy) datasets specified in the [tutorial](https://github.com/jvtalwar/GRIEVOUS/wiki/Data,-Preprocessing,-And-Formatting), go ahead and clone the GRIEVOUS github repository with `git clone https://github.com/jvtalwar/GRIEVOUS`. Note these datasets are located in the `Datasets/Sarcoidosis/` subdirectory. 

## Validating your GRIEVOUS' Tutorial Output:

For validation of a complete tutorial `grievous` run, we provide two files:

 1) **Tutorial_Intersecting_SNPs.tsv**:  A list of SNPs common to all 3 tutorial datasets, indexed according to **CHR**:**POS**:*sorted*(**REF**:**ALT**) notation.
    - *sorted*(**REF**:**ALT**) corresponds to alphabetical order of nucleotides for a position. For example 1:33:**A**:**C** as opposed to 1:33:**C:A**.
 2) **validate_tutorial.py**:  A script for comparing the equivalency of the output of your `grievous` tutorial run-through. Specifically this script takes in
                           the above Tutorial_Intersecting_SNPs.tsv file and the output of your `grievous intersect` run and assesses the equivalance in the outputs.
    - This script can also be used to assess the correctness of any given `grievous` run (post-intersect), provided a comprehensive set of intersecting SNPs in the 
      format of the Tutorial_Intersecting_SNPs.tsv file (this will need to be unique to your specific datasets of interest).


The run this validation script in accordance with the paths provided by the tutorial:

```bash
python validate_tutorial.py -i PATH_TO_GRIEVOUS_TUTORIAL_INTERSECTING_SNPs -g PATH_TO_GRIEVOUS_INTERSECT_OUTPUT
```

 - For example, according to the tutorial we would update the above PATHs as the following:
    - `PATH_TO_GRIEVOUS_TUTORIAL_INTERSECTING_SNPs=~/GRIEVOUS_Tutorial/Tutorial_Intersecting_SNPs.tsv`
    - `PATH_TO_GRIEVOUS_INTERSECT_OUTPUT=~/GRIEVOUS_Tutorial/IntersectingVariants/AllIntersectingSNPs.tsv`
 - These paths should be updated (as needed) according to any changes or modifications you made to your tutorial run through (in terms of paths and/or the output name of `grievous intersect`)
 - If you observe the below message, all went well with your `grievous` tutorial run through!
    ```
    PASS!!!
    grievous intersect output matches independently identified intersecting snps.
    ``` 