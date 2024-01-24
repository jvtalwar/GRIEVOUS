import grievous
import argparse
import logging


def main():
    logger = logging.getLogger()
    handler = logging.StreamHandler()

    parser = argparse.ArgumentParser()

    description = """
      _____ _____  _____ ________      ______  _    _  _____ 
     / ____|  __ \|_   _|  ____\ \    / / __ \| |  | |/ ____|
    | |  __| |__) | | | | |__   \ \  / / |  | | |  | | (___  
    | | |_ |  _  /  | | |  __|   \ \/ /| |  | | |  | |\___ \ 
    | |__| | | \ \ _| |_| |____   \  / | |__| | |__| |____) |
     \_____|_|  \_\_____|______|   \/   \____/ \____/|_____/

    Generalized Realignment of Innocuous and Essential Variants Otherwise Utilized as Skewed

    GRIEVOUS is a command-line tool for ensuring cross-cohort variant consistency. Specifically, 
    for any given cohort, GRIEVOUS will identify, orient, and realign all valid biallelic variants 
    from a chromosomal-level plink2 variant information file (i.e., .pvar) or summary statistic file to 
    a user-specific variant database. Upon database alignment completion for all datasets of interest, 
    GRIEVOUS returns a comprehensive, consistently oriented and indexed feature set (across all 
    datasets), which can be used for the user's task of choice (e.g., ML training, validation, and
    testing; PRS construction and validation, etc.). 
    
    Succinctly, GRIEVOUS organizes the variant galaxy to the standards of the trade federation.
    """

    parser = argparse.ArgumentParser(description = description, formatter_class = argparse.RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers(help = "Must be either: 'realign', 'merge', 'intersect', 'list_dbs', 'records', 'delete_database' (or of course 'general')", required = True, dest = "operation")

    realignParser = subparsers.add_parser("realign", help = "Generate realignment files for genomic (i.e., .pvar) or summary statistic (.ssf) files.", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    realignParser.add_argument("--file", "-f", type = str, required = True, help = "Path to file to be grievous aligned.")
    realignParser.add_argument("--database", "-d", type = str, default = None, help = """User defined alias for a grievous database/dictionary. If None (default) grievous
                            will check internal settings to see if a single registered database is established, and if so will employ that for realignment. If more than 
                            one database is registered, the user must specify a database alias. Database aliases can be accessed with 'grievous list_dbs'.""")
    realignParser.add_argument("--file_type", "-t", type = str, help = "File type in {'pvar', 'ssf'} when automatic detection of file extension fails (i.e., when file does not end with .pvar or .ssf).")
    realignParser.add_argument("--write_path", "-w", type = str, default = None, help = """Directory to which grievous realigned files will be written. 
                            If none provided, grievous will write aligned files to a new directory GrievousAlignedFiles in the directory where --file (-f) was provided.""")
    realignParser.add_argument("--mapping", "-m", type = str, default = None, help = "Path to a mapping file which registers columns of file(s) to be aligned to the grievous standard.")
    realignParser.add_argument("--comment_characters", "-c", nargs = '+', type = str, default = [], help = """All comment characters at the header of a file to be skipped. Accepts multiple
                            arguments in a space separated manner. Comment characters including a # should be wrapped in quotes (e.g., '##').""")
    realignParser.add_argument("--verbose", "-v", action = "store_true", default = False, help = "Flag that increases the output from the realignment process.")
    realignParser.add_argument("--return_all", "-r", action = "store_true", default = False, 
                            help = "Whether to return all columns of the original pvar/ssf file. Default: False.")
    realignParser.add_argument("--shutoff_db_update", "-s", action = "store_true", default = False, 
                            help = """Flag for whether database/dictionary updates should not be performed. To ensure grievous correctness this flag should only (if ever) be used on the final cohort/ssf 
                            for all datasets attempting to be realigned.""")


    mergeParser = subparsers.add_parser("merge", help = "Merge CHR level reports to composite reports. For ssf files, will also merge all CHR level grievous formatted ssfs to a composite ssf. Called after realign.", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    mergeParser.add_argument("--aligned_files", "-a", type = str, required = True, help = "Path to a completed grievous realign write path/out directory.")
    mergeParser.add_argument("--update_id", "-u", action = "store_true", default = False, help = """Flag for whether to update the IDs of a merged SSF with the grievous realigned IDs. Without this flag, merged SSF will
                            retain the old SSF IDs under the ID column, and the file will be indexed with the grievous realigned IDs. Flag has no effect on .pvar realigned directories.""")


    intersectParser = subparsers.add_parser("intersect", help = "Given n grievous realigned cohorts, generates an intersected feature set (i.e., all biallelic SNPs that are common to all n cohorts). Called after merge for each cohort.", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    intersectParser.add_argument("--write_path", "-w", type = str, required = True, help = "Directory to which grievous intersect output file will be written.")
    intersectParser.add_argument("--cohorts", "-c", nargs = '+', type = str, required = True, help = "All grievous realigned and merged directories desired for intersection. Accepts multiple arguments in a space separated manner.")
    intersectParser.add_argument("--exclude_biallelic_duplicates", "-xd", action = "store_true", default = False, help = """Flag as to whether to omit grievous identified biallelic SNPs that were identified as duplicated in any given cohort. 
                                By default these SNPs are included/returned by grievous intersect.""")
    intersectParser.add_argument("--exclude_palindromic", "-xp", action = "store_true", default = False, help = """Flag as to whether to omit grievous identified biallelic SNPs that are also palindromic SNPs.
                                By default these SNPs are included/returned by grievous intersect.""")
    intersectParser.add_argument("--out_name", "-n", type = str, default = "Grievous_Intersecting_SNPs", help = "Name of file returned by grievous intersect. Default: Grievous_Intersecting_SNPs")


    recordParser = subparsers.add_parser("records", help = """Identify all cohorts aligned and updated to a given grievous database/dictionary. Returns whole database records or a chromosomal subset
                                        of the database records when --chr parameter is employed.""", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    recordParser.add_argument("--database", "-d", type = str, default = None, help = """Alias for an already existing grievous database/dictionary. If None (default) grievous
                            will check to see if a single registered database is established, and if so it will report records for this database.""") 
    recordParser.add_argument("--chr", "-c", nargs = '+', type = str, default = None, help = "Report database records for specific chromosomes.") 
    recordParser.add_argument("--out", "-o", type = str, default = None, help = "Output file path to write database records. Default: None (i.e., write to current stream such as console).")


    listDictParser = subparsers.add_parser("list_dbs", help = "List all user-specified databases/dictionaries in your grievous library.", formatter_class = argparse.ArgumentDefaultsHelpFormatter)


    deleteDatabaseParser = subparsers.add_parser("delete_database", help = "Permanently delete an existing grievous database, including all variant orientations and records.", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    deleteDatabaseParser.add_argument("--database", "-d", type = str, required = True, help = "Alias for an already existing grievous database/dictionary.")


    generalParser = subparsers.add_parser("general", help = "Summon the general... if you dare.", formatter_class = argparse.ArgumentDefaultsHelpFormatter)


    args = parser.parse_args()

    if args.operation == "realign":
        if args.verbose:
            formatter = logging.Formatter("%(asctime)s - %(message)s")
            handler.setFormatter(formatter)
            logger.setLevel(logging.INFO)
            
        else:
            logger.setLevel(logging.WARNING)

        logger.addHandler(handler) 

        grievous.processes.Realign(file = args.file, dictionary = args.database, write_path = args.write_path, mapping = args.mapping, 
                                file_type = args.file_type, return_all = args.return_all, comment_characters = args.comment_characters, 
                                shutoff_dict_update = args.shutoff_db_update)
        
    elif args.operation == "merge":
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)

        grievous.processes.Merge(realignment_directory = args.aligned_files, updateIDs = args.update_id)

    elif args.operation == "intersect":
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)

        grievous.processes.Intersect(write_path = args.write_path, cohorts = args.cohorts, exclude_biallelic_duplicates = args.exclude_biallelic_duplicates,
                                    omit_palindromic = args.exclude_palindromic, file_output_name = args.out_name)

    elif args.operation == "records":
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)
        grievous.processes.utils.Records(dbAlias = args.database, chromosomes = args.chr, out = args.out)

    elif args.operation == "list_dbs":
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)
        grievous.processes.utils.ListDicts()

    elif args.operation == "delete_database":
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)
        grievous.processes.utils.DeleteDatabase(dbAlias = args.database)

    else:
        logger.setLevel(logging.INFO)
        logger.addHandler(handler) 
        grievous.processes.utils.GeneralGrievous()    

if __name__=="__main__":
	main()


