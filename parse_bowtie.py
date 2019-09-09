# By Robert W. Newberry, PhD
# Last Update: 190909 15:24
#
# Purpose: create a barcode:variant dictionary from Bowtie2 mappings
#
# Inputs
#    1) A Bowtie2 output file (SAM file)
#	    Obtained by aligning the output of merge.py against reference_sequences.fasta
#	    **The output file should be renamed 'bowtie_output.sam'
#    2) A list of barcodes that are well-represented in high-depth sequencing
#	    **Each barcode should be on a separate line
#	    **The list should be named 'barcode_list.txt'
#
# Outputs
#    1) A pickle dictionary relating barcodes to protein variants
#           Keys are 20bp barcodes
#           Values should be tuples, the first element being the position in sequence,
#               the second being the amino acid in one letter code
#               **WT is represented as (0,'WT')
#
# Environment: Python 2.7


import csv
from collections import Counter
import cPickle as pickle


# this function takes the output from running Bowtie on XXXX and returns the associated barcode:variant pairs for each read
def parse_bowtie():

# start by trimming the sam file to just the names of the query and reference sequences, which massively improves runtime
    with open('bowtie_output.sam', 'r') as fin, open('bowtie_output_trimmed.sam','w') as fout:
        sam_reader=csv.reader(fin,delimiter='\t',quoting=csv.QUOTE_NONE)
        sam_writer=csv.writer(fout,delimiter='\t',quoting=csv.QUOTE_NONE)
        for row in sam_reader:
            new_row=row[:5]
            sam_writer.writerow(new_row)

    with open('bowtie_output_trimmed.sam','r') as bowtie_file, open('bowtie_matches.txt','w') as fout: # open the variant and output files
        sam_reader=csv.reader(bowtie_file,delimiter='\t',quoting=csv.QUOTE_NONE) # create a csv reader for the variant sam file
        for row in sam_reader: # for each entry
            barcode=row[0][-20:] # the barcode is the last 20 characters of the query
            if row[2][-5] == row[2][-1]: # if the WT and mutant amino acids are the same,
                variant_call=(0,'WT') # then the barcode maps to WT
            else: # if not,
                variant_call=(int(row1[2][23:26]),row1[2][-1]) # then the variant is (position, AA)
            fout.write(barcode+','+str(variant_call)+'\n') # write the ouput to csv


# this function compiles all of the variants mapping to a particular barcode
def collapse():
    mappings={} # create a dictionary to house the variant mappings as mappings[barcode]=[variant1,variant2,etc.]
    with open('bowtie_matches.txt','r') as fin: # open the parse Bowtie results
        for line in fin: # for each barcode:variant pair
            try:
                mappings[line[:20]].append(line[21:-1]) # add the variant to the list of variants associated with that barcode
            except KeyError: # unless that barcode is not yet in the dictionary
                mappings[line[:20]]=[line[21:-1]] # in which case, create an entry in the dictionary for that barcode
    with open('bowtie_barcode_mappings.pkl','w') as fout:
        fout.write(pickle.dumps(mappings)) # write the ouput to a Pickle file


# this function identifies barcodes that reliably map to the same variant and creates a dictionary
def consensus():
    mappings=pickle.load(open('bowtie_barcode_mappings.pkl','r')) # open the dictionary of barcode:variant mappings
    best_maps={} # create a new dictionary for high-confidence barcode:variant mappings
    with open('barcode_list.txt', 'r') as barcode_list: # open the list of well-represented barcodes from high-depth sequencing
        for line in barcode_list: # for each well-represented barcode, 
            barcode=line[0:20]
            try:
                two_best_variants=Counter(mappings[barcode]).most_common(2) # find the two variants that most commonly map to a barcode
                if two_best_variants[0][1] > 1: # only consider mappings with at least two independent reads
            	    try:
                	if float(two_best_variants[0][1])/float(two_best_variants[1][1]) >= 2.0: # if the most common mapping is at least twice as common as the second most common mapping,
                    	    best_maps[barcode]=two_best_variants[0][0] # then take the most common mapping as high-confidence and write it to the new dictionary
              	    except IndexError: # if one a single variant maps to a barcode,
                        best_maps[barcode]=two_best_variants[0][0] # then add that barcode to the high-confidence dictionary
                except KeyError: # catch well-represented barcodes that weren't captured during long-read sequencing
                    donothing=0
    with open('consensus_dictionary.pkl','w') as fout: # write the high-confidence mappings to a new dictionary for use in barcode counting
        fout.write(pickle.dumps(best_maps))


parse_bowtie()
collapse()
consensus()
