# script to blast HIV DNA sequences to generate source sequences
# source represent sequences that were imported to the population being studied
# CHANGE DIRECTORY AND FILE NAME ACCORDINGLY - I STILL HAVE TO ADD A VARIABLE TO HANDLE THAT

from miscellaneous import Misc
import os
import csv



if __name__ == "__main__":

    # CHANGE FILENAME HERE
    # directory and filename containing the genbank output (after running blastn)
    gb_output = "/Your/Path/Here/SA_env_ali_output.xml"

    # CHANGE FILE NAME HERE
    # directory name to save new results
    dirname = "/Your/Path/Here/Blast_source/"

    # CHANGE NAME HERE
    # filename to save results
    filename = "env_ali_id_bestMatch.csv"

    misc = Misc()

    os.chdir(dirname)

    results = open(gb_output)

    seq_id_str = misc.check_blast_records1(results)
    results.close()
    print seq_id_str

    # Save results
    with open(filename, 'wb') as f:
        writer = csv.writer(f)
        for val in seq_id_str:
            writer.writerow([val])

