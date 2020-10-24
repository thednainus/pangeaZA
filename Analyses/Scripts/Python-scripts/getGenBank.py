# script to read a csv file containing one column for GenBank accession numbers
# and get the genbank file type for each of this accession numbers
# this information will be later used to have check country for each of the sequence to use as "source"
# in phylodynamic analysis

# You should always CHANGE the DIRECTORY NAME and FILENAMES below.
# You should also CHANGE country name on line 54


from miscellaneous import Misc
import csv
from Bio import SeqIO



if __name__ == "__main__":

    # CHANGE FILENAME HERE
    # directory and filename containing the list of GenBank accession numbers as csv file
    gb_ids = "/Your/Path/Here/bestMatchIDs/SA_env_ali_id_3bestMatches_genbank.csv"

    # CHANGE FILE NAME HERE
    # directory and filename to save new results
    gb_output_file = "/Your/Path/Here/Blast_source/South_Africa/GenBank_files/SA_src_env_ali_3bestMatches.gb"

    # directory and filename to save metadata for source sequences (accession number, country and year of HIV sample collection)
    src_metadata_filename = "/Your/Path/Here/Blast_source/South_Africa/Metadata/SA_src_env_ali_metadata.csv"

    # directory and filename to save fasta sequence from the selected genbank accesion numbers that was chose as best hits
    # not from South Africa
    fasta_output = "/Your/Path/Here/Blast_source/South_Africa/FASTA_files/SA_src_env_ali.fasta"

    misc = Misc()

    # open csv file and save elements as a list. However, this creates a list os lists
    with open(gb_ids, 'rb') as f:

        reader = csv.reader(f)
        genbankIDs_list = [r for r in reader]

    # in this chink of code I convert a list of list into a flat list to be used to get the genbank files
    # for the GenBank accession numbers
    gb_flat_list = []
    for sublist in genbankIDs_list:
        for item in sublist:
            gb_flat_list.append(item)

    misc.download_genbank(gb_flat_list, gb_output_file)
    gb_records = SeqIO.parse(open(gb_output_file, "r"), "genbank")
    #CHANGE COUNTRY NAME HERE
    #This line is teoretically not necessary anymore because now I have a 
    # specific database in which I exclude sequences from countries that I am 
    # interested. However, I did not modify it in the script, and this line can
    # also be used to make sure I don't have sequences not outside from country 
    # of interest in my FASTA sequence file.
    country_name = "south africa"
    source = misc.get_source_genbank(gb_records, country_name, fasta_output, src_metadata_filename)


