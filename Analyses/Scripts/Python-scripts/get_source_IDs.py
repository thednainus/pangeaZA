# Script to get best matches (PANGEA or GenBank IDs) to generate source sequences using script getGenBank.py
# (for GenBank sequences) or R script for PANGEA sequences.
# source represents sequences that were imported to the population being studied

# CHANGE DIRECTORY AND FILE NAME ACCORDINGLY - I STILL HAVE TO ADD A VARIABLE TO HANDLE THAT

from miscellaneous import Misc
import os
import csv



if __name__ == "__main__":

    # CHANGE FILENAME HERE
    # directory and filename containing the genbank output (after running blastn)
    gb_output = "/Your/Path/Here/SA_env_ali.xml"

    # CHANGE FILE NAME HERE
    # directory name to save new results
    dirname = "/Your/Path/Here/bestMatchIDs/"

    # CHANCE NAME HERE
    # filename = only best ID matches (it can be PANGEA or GenBank ids)
    filename = "SA_env_ali_id_3bestMatches.csv"
    # filename2 = the same as filename, but results here are in triplets: 
    # query sequence, best hit, order of best hit
    filename2 = "SA_env_ali_id_3bestMatches_triplets.csv"

    misc = Misc()

    os.chdir(dirname)

    results = open(gb_output)

    seq_id_str = misc.check_blast_records3(results, n=3)
    results.close()
    print seq_id_str

    # Save results (non-identical blast matches for all query sequences
    with open(filename, 'wb') as f:
        writer = csv.writer(f)
        for val in seq_id_str[0]:
            writer.writerow([val])


    # Save results. This is the same as above csv.write. The difference is that 
    #I am saving triplets of query sequence, and best hit. This might be helpful 
    # if by any chance one decides do remove query sequences from the alignment, 
    # in which case this list can be used to remove the correspondent 
    # CGR sequence if exclusive to a specifi CGR sequence.
    with open(filename2, 'wb') as f2:
        writer = csv.writer(f2, delimiter=",")
        for val2 in seq_id_str[1]:
            writer.writerow(val2)

