import os 
import logging
import random

logger = logging.getLogger(__name__)

'''
About: Method to check if .grievous directory has been created, which serves as the grievous 
dictionary/database library. If it does not exist this method creates it. Returns the path
to the .grievous directory.

Input(s): None
Output(s): grievousDir: String corresponding to the .grievous directory which contains all databases.
'''
def _QueryGrievousInstall():
    homeDir = os.path.expanduser('~')
    grievousDir = os.path.join(homeDir, ".grievous")

    if not os.path.isdir(grievousDir):
        os.makedirs(grievousDir)
        assert os.path.isdir(grievousDir), f"Unable to create .grievous folder in {homeDir}. Exiting..."

    return grievousDir

'''
About: Method to display all user-created grievous databases.

Input(s): None
Output(s): None
'''
def ListDicts():
    grievousDir = _QueryGrievousInstall()

    #get database/dictionary aliases
    grievousLibrary = [el for el in os.listdir(grievousDir) if os.path.isdir(os.path.join(grievousDir, el))]

    logger.info("\nCurrent grievous database aliases:\n" + "-"*37)
    for dictionary in grievousLibrary:
        logger.info(dictionary)

    logger.info("")

    return None 

'''
About: Method to identify all palindromic SNPs.

Input(s): candidateSNPs: Iterable (i.e., set, list, tuple) of SNPs
Output(s): palindromicSNPs: set of all palindromic SNPs in candidate SNPs
'''
def _IdentifyPalindromic(candidateSNPs):
    palindromicPairs = {"A:T", "C:G"}
    palindromicSNPs = set()
    for snp in candidateSNPs:
        snpPair = ":".join(sorted(snp.split(":")[-2:]))
        if snpPair in palindromicPairs:
            palindromicSNPs.add(snp)

    return palindromicSNPs

'''
About: This method generates a record report for a given database. Specifically
it returns all cohorts/files that have been grievous realigned to a given database/dictionary
and employed for dictionary/database updates. This method returns records to logger 
(i.e., user stream) in the absence of out (i.e., out is None).

Input(s): 1) dbAlias: String corresponding to the grievous dictionary/database for which the user is attempting to access records.
          2) chromosomes: List of strings corresponding to a subset of particular chromosomes for which to access records (as opposed to all).
          3) out: String corresponding to a path (including file name and extension) to which to write files.

Output(s): None
Write(s):  recordReport - A report of all aligned and updated files aligned to a given database. If chromosomes is provided 
           then report is limited to those chromosomes of interest. Otherwise returns records for all chromosomes.
'''
def Records(dbAlias, chromosomes, out):
    if out: #Validate out path if user provided one:
        if os.path.dirname(out) != '': 
            assert os.path.isdir(os.path.dirname(out)), f"--out parent directory {os.path.abspath(os.path.dirname(out))} does not exist. Provide valid path to --out."

        if os.path.isdir(out): #user provided a directory path with no file name
            out = os.path.join(out, "GRIEVOUS_RECORDS.txt")

    grievousDir = _QueryGrievousInstall()

    #get database/dictionary aliases
    grievousLibrary = [el for el in os.listdir(grievousDir) if os.path.isdir(os.path.join(grievousDir, el))]

    #check database alias or if none provided then ensure only single db exists
    grievousRecordsDir = ""

    if dbAlias not in grievousLibrary:
        if dbAlias is None:
            assert len(grievousLibrary) == 1, """--database not provided and either more than one database exists or none have been created. 
            Provide a current database/dictionary alias to access corresponding records. Database aliases can be found with 'grievous list_dbs'."""

            grievousRecordsDir = os.path.join(grievousDir, grievousLibrary[0], "Records")
            dbAlias = grievousLibrary[0]
        
        else: #database alias provided by user does not exist
            raise Exception(f"--database {dbAlias} not found among grievous databases. Current grievous databases can be found with 'grievous list_dbs'.")

    else:
        grievousRecordsDir = os.path.join(grievousDir, dbAlias, "Records")

    chromosomalIteration = [str(i) for i in range(1,23)] + ["X", "Y", "MT"]
    
    if chromosomes: # user provided chromosomal subset
        #validate to ensure is a subset of chromosomalIteration
        for chromosome in chromosomes:
            assert chromosome in chromosomalIteration, f"--chr {chromosome} unknown. Valid --chr inputs are 1-22, X, Y, MT."

        chromosomalIteration = chromosomes

    #Add line at top of report indicating the database records reported
    if out:
        with open(out, "a") as recordReport:
            recordReport.write(f"Reporting Records For Database:  {dbAlias} \n\n")

    else:
        logging.info(f"Reporting Records For Database:  {dbAlias} \n")

    for chromosome in chromosomalIteration:
        with open(os.path.join(grievousRecordsDir, f"CHR_{chromosome}.aligned.record.txt"), "r") as alignedRecord:
            alignedContent = alignedRecord.readlines()

        with open(os.path.join(grievousRecordsDir, f"CHR_{chromosome}.updated.record.txt"), "r") as updatedRecord:
            updatedContent = updatedRecord.readlines()
        
        if out: #write to file if user provided one
            with open(out, "a") as recordReport:
                recordReport.write(f"CHR {chromosome} Records:\n\nALIGNED\n\n")
                for line in alignedContent:
                    recordReport.write(line)

                recordReport.write("\nUPDATED DATABASE:\n\n")

                for line in updatedContent:
                    recordReport.write(line)

                recordReport.write("\n")

        else:
            logger.info(f"CHR {chromosome} Records:\n\nALIGNED\n")
            for line in alignedContent:
                logger.info(line)

            logger.info("UPDATED DATABASE:\n") 

            for line in updatedContent:
                logger.info(line)
            

    return None

'''
About: This method takes in a database alias and if it exists deletes the database from the 
grievous library.

Input(s): dbAlias: String corresponding to the grievous dictionary/database the user is attempting to delete.
Output(s): None
'''
def DeleteDatabase(dbAlias):
    grievousDir = _QueryGrievousInstall()

    #get database/dictionary aliases
    grievousLibrary = [el for el in os.listdir(grievousDir) if os.path.isdir(os.path.join(grievousDir, el))]

    if dbAlias not in grievousLibrary:
        raise Exception(f"--database {dbAlias} not found among grievous databases. Current grievous databases can be found with 'grievous list_dbs'.")
    
    databaseToBeDeleted = os.path.join(grievousDir, dbAlias)
    recordsToRemove = os.path.join(databaseToBeDeleted, "Records")

    #delete records
    if os.path.isdir(recordsToRemove): #handle case of user interruption and entire Records is deleted, but not the entire database
        for record in os.listdir(recordsToRemove):
            toRemove = os.path.join(recordsToRemove, record)
            os.remove(toRemove)

        os.rmdir(recordsToRemove)

    #delete all chromosomal level orientation databases
    for chromDatabase in os.listdir(databaseToBeDeleted):
        toRemove = os.path.join(databaseToBeDeleted, chromDatabase)
        os.remove(toRemove)
    
    os.rmdir(databaseToBeDeleted)

    return None

'''
About: Every story needs a good old-fashioned villain... 
(Translation: What good is naming a tool grievous if you can't
have fun with it?). Queue Revenge of the SNP!

Input(s): None
Output(s): None
'''
def GeneralGrievous():
    quotes = ["If a variant does not appear in our records, it does not exist... yet.",
              "Your variants will make a fine addition to my collection!",
              "I will deal with this variant slime myself...",
              "I don't like SNPs. They're coarse and rough and irritating and they get everywhere.",
              "Wait a minute. How did this happen? We're smarter than this.\nApparently not.",
              "General Grievous, you're shorter than I expected.",
              "Hello there!",
              "Army or not, you must realize you are doomed!",
              "Only a SNP deals in absolutes. I will do what I must.",
              "SNP lords are our speciality.",
              "No giving up general jar jar. Weesa think of something.\nMy give up. My give up.",
              "Activate ray shields!",
              "Don't bother with them. Keep the SNPs in orbit!"]
    
    chosenOne = random.choice(quotes)

    grievousBanner = ["@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@","@@@@@@@@@@@@@@@@@@::..::@@@@@@@@@@@@@@@@@@",
                      "@@@@@@@@@@@@@@@@@:......::@@@@@@@@@@@@@@@@","@@@@@@@@@@@- @@@:.:......::@@ @@@@@@@@@@@@",
                      "@@@@@@@@@@@@ @@@:.=:....-::@@:@@@@@@@@@@@@","@@@@@@@@@@@@.@@-::=.....--:@@.@@@@@@@@@@@@",
                      "@@@@@@@@@@@. @@-.:-.....:=:=@  @@@@@@@@@@@","@@@@@@@@@@ ..@+-..=.....:::-@.  @@@@@@@@@@",
                      "@@@@@@@@@ ..:-#-..-.....:=--@..  @@@@@@@@@","@@@@@@@@ ..::--:.::.....:#::-::.  @@@@@@@@",
                      "@@@@@@@....::-=:.-:.:.:::#::+....  @@@@@@@","@@@@@@@...:*:-::.*..:...:#::--%:.. @@@@@@@",
                      "@@@@@@@:...@*-:::#::.::::#::::*#.  @@@@@@@","@@@@@@@@..+##-:::*:::::::#::--@@.. @@@@@@@",
                      "@@@@@@@@.:@#*:%@%@%:.::-*@%@%=@@:. @@@@@@@","@@@@@@@@:.-**:%@*+#*:::+@#*@%:#@..:@@@@@@@",
                      "@@@@@@@@:.:*#-%@@@#:::::*@@@%-*@..@@@@@@@@","@@@@@@@@..:=#==%@*::::::-=@%==*.. @@@@@@@@",
                      "@@@@@@@@=.::%#-----::::---=-:%+:. @@@@@@@@","@@@@@@@@@-:.-%=-==-:::----==-+:. .@@@@@@@@",
                      "@@@@@@@@@@-::=*++=-::::-===++-:::@@@@@@@@@","@@@@@@@@@@@::-+@%=-::::--=##=-::@@@@@@@@@@",
                      "@@@@@@@@@@@+::=%@#-::::--*@%-::@@@@@@@@@@@","@@@@@@@@@@@@-:-%@#-::::--*@#:::@@@@@@@@@@@",
                      "@@@@@@@@@@@@@::#@@-::::--%@+::@@@@@@@@@@@@","@@@@@@@@@@@@++::%@-:::--=%%::-@@@@@@@@@@@@",
                      "@@@@@@@@@@@+**:.%@-:----=@%.:%@@@@@@@@@@@@","@@@@@@@@@@*%*@=.@%-=*+===#@.-%%@@@@@@@@@@@",
                      "@@@@@@@@@++%%@@.#@-=*+==-@@.%@#+@@@@@@@@@@","@@@@@@@@-=**@@@-#@:+*++=-@%.@@##-@@@@@@@@@",
                      "@@@@@@@@==#%@@@%*%++***==##-@@@%=-@@@@@@@@","@@@@@@@==+%#@@@@:*:***++-*#@@@%%+=-@@@@@@@",
                      "@@@@@@@==#@%@@@@#=:++#+*-*@@@@#@%+-@@@@@@@","@@@@@@-=*%@@*@@@%=:=****-*@@@@%@@*-:@@@@@@",
                      "@@@@+*-=*@@@#@@%@+:*####-#@@@%@@@%-:*@@@@@"]

    logger.info("\n")
    for line in grievousBanner:
        logger.info(line)

    logger.info(f"\n\n{chosenOne}\n")

    return None