
#@author: James V. Talwar

from ..GenomicObjects import Pvar, SSF
import os
import pandas as pd
import logging
import time
from .utils import _QueryGrievousInstall

logger = logging.getLogger(__name__)

'''
About: Helper method to create a new grievous database/dictionary directory.

Input(s):  1) newDictPath: String corresponding to the directory where empty database/dictionary records will be created to be updated.
           2) dictStorage: String corresponding to database storage method. Default: 'parquet'
Output(s): None
'''
def _Create_Dictionary(newDictPath, dictStorage = "parquet"):
    if os.path.isdir(newDictPath): #chromosome-level paralleleization - if database is in process of (or has finished) being created don't need to re-create/overwrite db
        return None

    #Create new dictionary folder
    os.makedirs(newDictPath, exist_ok = True)
    os.makedirs(os.path.join(newDictPath, "Records"), exist_ok = True)

    #Create empty record files and empty dictionary files 
    chromosomes = [str(i) for i in range(1,23)] + ["X", "Y", "MT"] #Grievous standard for chromosomes

    for chromosome in chromosomes:
        # all the files aligned to given database 
        with open(os.path.join(newDictPath, f"Records/CHR_{chromosome}.aligned.record.txt"), "a") as alignedRecord: 
            pass

        # all the files which were used to update the given database. If --shutoff_dict_update is never used then is equivalent to aligned
        with open(os.path.join(newDictPath, f"Records/CHR_{chromosome}.updated.record.txt"), "a") as updatedRecord: 
            pass

        if dictStorage == "parquet":
            parquetWrite = os.path.join(newDictPath, f"CHR_{chromosome}.parquet")
            if not os.path.exists(parquetWrite):
                pd.DataFrame([], columns = ["CHR", "POS", "REF", "ALT", "ID"]).to_parquet(parquetWrite)
            
        else:
            tsvWrite = os.path.join(newDictPath, f"CHR_{chromosome}.tsv")
            if not os.path.exists(tsvWrite):
                pd.DataFrame([], columns = ["CHR", "POS", "REF", "ALT", "ID"]).to_csv(tsvWrite, sep = "\t")
            

    return None

'''
Input(s): dbAlias: String corresponding to an already created or novel dictionary/database alias. 
Outputs: grievousDictDir: String corresponding to the full path of the desired grievous dictionary/database directory.
'''
def _Check_Dictionary(dbAlias):
    grievousDir = _QueryGrievousInstall()

    #check among the dictionary aliases
    grievousDictAliases = [el for el in os.listdir(grievousDir) if os.path.isdir(os.path.join(grievousDir, el))]

    if dbAlias not in grievousDictAliases:
        if dbAlias is None: #Handle automatically when only one dictionary standard exists.
            assert len(grievousDictAliases) == 1, """--database not provided and either more than one database exists or none have been created. 
            Provide a current database/dictionary alias or provide a new one to create a new empty database. Database aliases can be found with grievous list_dbs."""

            grievousDictDir = os.path.join(grievousDir, grievousDictAliases[0])

        else:
            logger.warning(f"CREATING A NEW GRIEVOUS DICTIONARY WITH ALIAS: {dbAlias}")
            
            grievousDictDir = os.path.join(grievousDir, dbAlias)
            _Create_Dictionary(grievousDictDir)

    else:         
        grievousDictDir = os.path.join(grievousDir, dbAlias)

    # Ensure all chromosomal dictionaries/dbs are written; instances of creation and file non-existence can occur when initializing a DB and running in parallel. 
    # To prevent a given chromosome from attempting to access a db before creation (while in progress) sleep to allow for all chromosome creation.
    if not (set([el.split(".")[0].split("_")[1] for el in os.listdir(grievousDictDir) if "CHR" in el]) == set([str(i) for i in range(1,23)] + ["X", "Y", "MT"])):
        time.sleep(13)


    return grievousDictDir


'''
About: This method performs grievous realignment, by creating a PVAR/SSF object and then calling the instance methods. For full details of these methods see
package GenomicObjects.

Input(s): 1) file: String corresponding to the path of the file which the user is attempting to realign.
          2) dictionary: String corresponding to the dictionary/database to which the user is attempting to realign. 
          3) write_path: String corresponding to the directory where the user wants to write realignment files.
          4) mapping: String corresponding to the file path registering all mandatory non-grievous standard column names to the grievous standard.
          5) file_type: String in {'pvar', 'ssf'} corresponding to the desired genomic file type to be realigned.
          6) return_all: Boolean corresponding to whether to return all columns of the underlying file, or the grievous subset.
          7) comment_characters: List of all comment characters corresponding to lines to skip at the beginning of a genomic file.
          8) shutoff_dict_update: Boolean corresponding to whether to skip the dictionary update step. If there is ever any ambiguity 
          as to whether to use this or not (or whether one will align future cohorts), then this parameter SHOULD ALWAYS BE FALSE!

Output(s): None
Write(s): Grievous PVAR/SSF output files (for full details see Pvar.Write() or SSF.Write()).
'''
def Realign(file, dictionary, write_path, mapping, file_type, return_all, comment_characters, shutoff_dict_update):
    assert os.path.exists(file), f"File path {file} does not exist. Exiting..."
    
    if write_path is not None: #If user does not pass in write_path will attempt to create a write directory: GrievousAlignedFiles
        assert os.path.isdir(write_path), f"--write_path {write_path} does not exist as a directory. Please retry with a valid write_path. Exiting..."
    
    #load the file and instantiate proper grievous genomic object
    if file_type not in {"pvar", "ssf"}:
        logger.warning("Provided file type not in {'pvar', 'ssf'}. Attempting automatic detection of file type...")
        
        extension = file.split(".")[-1].lower()
        if extension not in {"pvar", "ssf"}:
            logger.critical("Automatic detection of file type failed. Extension does not end with pvar or ssf.")
            raise Exception("Unable to resolve file type. Please provide --file_type. Options include pvar and ssf. Exiting...")
        
        else:
            logger.warning(f"Automatic detection of file type resolved to {extension}\n")
            file_type = extension    

    if file_type == "pvar":
        gFile = Pvar(filePath = file, mappingPath = mapping, commentCharsToSkip = comment_characters)
    
    else:
        gFile = SSF(filePath = file, mappingPath = mapping, commentCharsToSkip = comment_characters)

    #Identify/Make necessary grievous database/dictionary directory and load chromosome-level database/dictionary
    grievousDictionaryPath = _Check_Dictionary(dictionary)
    dictionaryExtension = ""

    try:
        chromosomeDictionary = pd.read_parquet(os.path.join(grievousDictionaryPath, f"CHR_{gFile.Chrom}.parquet"))
        chromosomeDictionary["CHR"] = chromosomeDictionary["CHR"].astype(str)

        dictionaryExtension = "parquet"

    except:
        logger.info("parquet chromosome database not found... Trying tsv")
        assert f"CHR_{gFile.Chrom}.tsv" in os.listdir(grievousDictionaryPath), f"No valid chromosomal database file found in grievous database {dictionary} for chromosome {gFile.Chrom}. Exiting..."
        chromosomeDictionary = pd.read_csv(os.path.join(grievousDictionaryPath, f"CHR_{gFile.Chrom}.tsv"), sep = "\t", index_col = 0, dtype={'CHR': str})

        dictionaryExtension = "tsv"


    #if user fails to provide write_path attempt to make a write directory for them in the folder where the file to be aligned is:
    if write_path is None:
        fileParentDir = os.path.dirname(file)
        write_path = os.path.join(fileParentDir, "GrievousAlignedFiles")
        logger.warning(f"No write_path provided. Writing aligned files to a new directory {write_path}")
        os.makedirs(write_path, exist_ok = True) 

    else:
        if (os.path.exists(os.path.join(write_path, "GRIEVOUS_Formatted", f"FormattedAndAlignedPvar_CHR{gFile.Chrom}.tsv"))) or (os.path.exists(os.path.join(write_path, "GRIEVOUS_Formatted", f"FormattedAndAlignedSSF_CHR{gFile.Chrom}.tsv"))):
            logger.warning(f"Previous GRIEVOUS_Formatted file detected in write_path {write_path}. If this directory was used for a previous realignment, current files will be overwritten.") 

    # Make reports, reorientation, and greivousFormatted folders for output files
    os.makedirs(os.path.join(write_path, "Reports"), exist_ok = True)
    os.makedirs(os.path.join(write_path, "GRIEVOUS_Formatted"), exist_ok = True)

    if file_type == "pvar":       
        os.makedirs(os.path.join(write_path, "Reorientation"), exist_ok = True)



    #Implement Realignment

    #Employ Clean:
    gFile.Clean()

    #Employ Orient:
    gFile.Orient(chromosomeDictionary)

    #Employ Align:
    gFile.Align(write_path)

    if not shutoff_dict_update: 
        #Update grievous database/dictionary     
        gFile.UpdateDict(chromosomeDictionary,  grievousDictionaryPath, dictionaryExtension)

        #Update .updated.record file
        with open(os.path.join(grievousDictionaryPath, f"Records/CHR_{gFile.Chrom}.updated.record.txt"), "a") as updateRecord:
            updateRecord.write(file + "\n")
        

    #Write realigned files:
    gFile.Write(write_path, return_all)

    #update .aligned.record file for the given grievous database
    with open(os.path.join(grievousDictionaryPath, f"Records/CHR_{gFile.Chrom}.aligned.record.txt"), "a") as alignedRecord:
        alignedRecord.write(file + "\n")
    
    return None 