from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO
import csv


class Misc(object):

    def create_dict_from_variables(self, key, value):

        variable_dict = {key: value}
        return variable_dict

    # get the first best match
    def check_blast_records1(self, results):

        count = 0

        # create list with IDs for pangea ID or GenBank accession numbers
        id_list = []

        blast_records = NCBIXML.parse(results)

        E_VALUE_THRESH = 1e-3
        for blast_record in blast_records:

            query = str(blast_record.query)
            query = query.split(".")
            print(query[0])
            print ('query name', blast_record.query)
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:

                    if hsp.expect < E_VALUE_THRESH:

                        if count < 1:

                            print('****Alignment****')
                            print('sequence:', str(alignment.hit_id))
                            seq = str(alignment.hit_id)
                            print(seq)
                            #if pangea, sequence id will be something on the format PG16-BW002502
                            #if sequence deposited in genbank it will be something in the format gb|AY713415|
                            seq_list = seq.split("|")

                            # so if it is a sequence deposited in genbank I am interested in the genbank number
                            if len(seq_list) > 1:
                                seq = seq_list[1]

                            if seq != query[0]:

                                #print("Genbank: ", GenBank[0])
                                print("GI: ", seq[1])

                                update_dict = self.create_dict_from_variables(query[0], seq)

                                dict.update(update_dict)

                                if seq not in id_list:
                                    id_list.append(seq)

                                count = count + 1
                                print (count)

            count = 0  # reset count to zero for the next blast query


        return id_list

    # Function to get best matches to check later, using script getGenBank.py, whether they are from a country different
    # from the country of interest, or whether country is not available
    # note that n will determine if you would like the best 3 matches or all 25 best matches (in this case
    # I did a blastn and saved best 25 matches. If one decide to save more matches, one can set the number to n to
    # the value of interest
    def check_blast_records3(self, results, n):
        count = 0

        # create list with IDs for pangea ID or GenBank accession numbers
        id_list = []

        # list to create pairs of sequences
        # First element will be Pangea ID query sequences and
        # second element will be Pangea ID or GenBank accession number for best hit
        pairs_list = []

        blast_records = NCBIXML.parse(results)

        E_VALUE_THRESH = 1e-3
        for blast_record in blast_records:

            query = str(blast_record.query)
            query = query.split(".")
            # print(query[0])
            # print ('query name', blast_record.query)
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= E_VALUE_THRESH:

                        # This is the first blast hit (first best match when count < 0 or count == 0)
                        if count < 1:

                            # print('****Alignment****')
                            # print('sequence:', str(alignment.hit_id))
                            seq = str(alignment.hit_id)
                            # print(seq)
                            # if pangea, sequence id will be something on the format PG16-BW002502
                            # if sequence deposited in genbank it will be something in the format gb|AY713415|
                            seq_list = seq.split("|")

                            # so if it is a sequence deposited in genbank I am interested in the genbank number
                            if len(seq_list) > 1:
                                seq = seq_list[1]

                            # check if first best match is identical to query name
                            # if it is, it will move to the next best hit
                            if seq != query[0]:

                                if seq not in id_list:
                                    id_list.append(seq)
                                    pairs = [query[0], seq, count+1]
                                    pairs_list.append(pairs)

                                # if seq already in list I should not append the same element to pairs_list
                                # but I should increase count anyways
                                count = count + 1
                                # print (count)

                        elif count < n and count != 0:

                            # print('****Alignment****')
                            # print('sequence:', str(alignment.hit_id))
                            seq = str(alignment.hit_id)
                            # print(seq)
                            # if pangea, sequence id will be something on the format PG16-BW002502
                            # if sequence deposited in genbank it will be something in the format gb|AY713415|
                            seq_list = seq.split("|")

                            # so if it is a sequence deposited in genbank I am interested in the genbank number
                            if len(seq_list) > 1:
                                seq = seq_list[1]

                            # check if first best match is identical to query name
                            # if it is, it will move to the next best hit
                            if seq != query[0]:

                                if seq not in id_list:
                                    id_list.append(seq)
                                    pairs = [query[0], seq, count+1]
                                    pairs_list.append(pairs)

                                # if seq already in list I should not append the same element to pairs_list
                                # but I should increase count anyways
                                count = count + 1

            count = 0  # reset count to zero for the next blast query


        return [id_list, pairs_list]


    def download_genbank(self, genbank_ids, gb_output_file):
        Entrez.email = "f.nascimento@imperial.ac.uk"
        handle = Entrez.efetch(db="nuccore", id=genbank_ids, rettype="gb", retmode="text")
        records = SeqIO.parse(handle, "gb")

        # save GenBak records to file
        output_file = open(gb_output_file, 'w')
        SeqIO.write(records, output_file, "genbank")




    # get gi for entries in which country is NOT as in "country_name"
    # it will also save as fasta the genbank entries for the source sequences
    def get_source_genbankCORRECT(self, gb_records, country_name, fasta_output, src_output):

        key = "country"

        # list of GenBank accession numbers that does not match "Senegal" as country
        # This will be the list containing the source GenBank accession numbers
        accessions_list = []
        country_list = []
        year_list = []

        count1 = 0
        count2 = 0

        output_file = open(fasta_output, "a")

        for record in gb_records:
            count1 += 1

            for f in record.features:
                count2 += 1

                if f.qualifiers.get(key) != None:

                    country = f.qualifiers.get(key)
                    print("country: ", country[0])
                    country = country[0].lower()

                    # sometimes country can come on the form of 'South Africa: Durban'
                    # split will allow to get just the country name
                    country_split = country.split(":")

                    if len(country_split) > 1:
                        country = country_split[0]

                    if country_name not in country:
                        # print(record.name)
                        accessions_list.append(record.name)
                        country_list.append(country)
                        if f.qualifiers.get("collection_date") != None:
                            year_list.append(f.qualifiers.get("collection_date")[0])
                        else:
                            year_list.append("N/A")

                        SeqIO.write(record, output_file, "fasta")

        output_file.close()
        self.save_source_metadata(src_output, accessions_list, country_list, year_list)

        return accessions_list




    # filename should include filename to save the metadata and path to save it
    # as for example: "/Users/user/Desktop/SOURCE_TESTE.csv"
    def save_source_metadata(self, filename, genbank_list, country_list, year_list):

        with open(filename, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            # write header to csv file
            writer.writerow(("GenBank", "Country", "Collection_year"))

            for genbank, country, year in zip(genbank_list, country_list, year_list):
                writer.writerow([genbank, country, year])


    # get gi for entries in which country is NOT as in "country_name"
    # it will also save as fasta the genbank entries for the source sequences
    def get_source_genbank(self, gb_records, country_name, fasta_output, src_output):
        #key = "country"

        # list of GenBank accession numbers that does not match "Senegal" as country
        # This will be the list containing the source GenBank accession numbers
        accessions_list = []
        country_list = []
        year_list = []

        count = 0

        output_file = open(fasta_output, "a")

        for record in gb_records:


            for f in record.features:

                # get just the first qualifier. It appears that only the first qualifier has information on country
                # usually a record will have more than 1 qualifier, and if I try to get country in a qualifier that
                # does not have information on country it will return None.
                if count == 0:

                    count += 1
                    which_country = f.qualifiers.get("country")
                    date = f.qualifiers.get("collection_date")

                    if which_country != None:
                        country = which_country[0].lower()
                    elif which_country == None:
                        country = "NA"

                    # sometimes country can come on the form of 'South Africa: Durban'
                    # split will allow to get just the country name
                    country_split = country.split(":")

                    if len(country_split) > 1:
                        country = country_split[0]

                    if country_name not in country:
                        # print(record.name)
                        accessions_list.append(record.name)
                        country_list.append(country)
                        if date != None:
                            year_list.append(date[0])
                        else:
                            year_list.append("NA")

                        SeqIO.write(record, output_file, "fasta")

            count = 0

        output_file.close()
        self.save_source_metadata(src_output, accessions_list, country_list, year_list)

        return accessions_list
