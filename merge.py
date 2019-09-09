# By Robert W. Newberry, PhD
# Last Update: 190909 14:51
#
# Purpose: merge long-read sequencing data for a-synuclein variants
#
# Inputs
#    1) FASTQ files for R1, R2, and I1, named r1.fastq, r2.fastq, i1.fastq, respectively
#       In this experiment, R1 is 300 bp, R2 is 130 bp, and I1 is 20 bp
#	Other designs can be accommodated by changing the read lengths
#
# Outputs
#    1) A FASTQ file containing merged DNA sequences
#	    The name contains both the cluster ID and the barcode sequence
#
# Environment: Python 2.7
#
# Known Issues
#     This script is specific to the design of a specific sequencing experiment;
#         other designs can be accommodate by changing the WT sequence and read lengtsh


import cPickle as pickle
from itertools import islice	


def calc_quality(phred_seq):
    overall_prob=1
    for letter in phred_seq:
        prob_correct=1-(10**(-(float((ord(letter)-33))/10.)))
        overall_prob*=prob_correct
    return overall_prob


def reverse_complement(seq):
    complement = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    seq = list(seq)
    seq.reverse()
    seqrc = ''
    for s in seq:
        seqrc = seqrc+complement[s]
    return seqrc


def reverse_complement_file(): # creates a fastq file for the reverse complement of each read
    counter=2
    with open('r2.fastq','r') as forward_file, open(filename[:-6]+'-rc.fastq','w') as reverse_file:
        for line in forward_file:
	    counter+=1
	    if (counter/4.).is_integer():
                reverse_file.write(reverse_complement(line[0:160])+'\n')
	    else:
		reverse_file.write(line)


def merge(): # merge the two reads and reconcile overlapping bases
    WT_sequence='ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGACAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCC'
    with open('r1.fastq','r') as read1_file, open('r2-rc.fastq','r') as read2_file, open('i1.fastq','r') as index_file, open('merged.fastq','w') as output_file:
	codes=islice(index_file,0,None,4)
	r1_reads=islice(read1_file,1,None,2)
	r1_scores=islice(read1_file,1,None,2)
        r2_reads=islice(read2_file,1,None,2)
        r2_scores=islice(read2_file,1,None,2)
	for line in codes:
	    output_file.write(line)
	    r1_read=r1_reads.next()
	    r1_score=r1_scores.next()
            r2_read=r2_reads.next()
            r2_score=r2_scores.next()
	    merged_read_list=list(r1_read[0:290]) # the merged read starts with the first 290 bases of R1
	    merged_score_list=list(r1_score[0:290]) # and the associated quality scores
	    for x in range(10): # for the region in which R1 and R2 overlap,
		if r1_read[x+290]==r2_read[x]: # if the reads agree, use that call
		    merged_read_list.append(r1_read[x+290])
		    merged_score_list.append(chr(max(ord(r1_score[x+290]),ord(r2_score[x]))))
		elif r1_read[x+290]==WT_sequence[x+290]: # if they disagree, but one indicates WT
		    merged_read_list.append(r1_read[x+290]) # call the base WT
		    merged_score_list.append(r1_score[x+290])
		elif r2_read[x]==WT_sequence[x]:
		    merged_read_list.append(r2_read[x])
		    merged_score_list.append(r2_score[x])
		else: # if they disagree and neither is WT, take the base with the higher quality score
		    if ord(r1_score[x+290]) > ord(r2_score[x]):
                        merged_read_list.append(r1_read[x+290])
                        merged_score_list.append(r1_score[x+290])
                    elif ord(r1_score[x+290]) < ord(r2_score[x]):
                        merged_read_list.append(r2_read[x])
                        merged_score_list.append(r2_score[x])
		    else:
			merged_read_list.append('N')
			merged_score_list.append(r2_score[x])
	    merged_read_list.append(r2_read[10:130])
	    merged_score_list.append(r2_score[10:130])
	    merged_read="".join(merged_read_list)
	    merged_score="".join(merged_score_list)
	    output_file.write(merged_read+'\n')
	    output_file.write('+\n')
	    output_file.write(merged_score+'\n')


def pair(): # takes the merged sequence file and the barcode file
    with open('merged.fastq','r') as allele_file, open('i1.fastq','r') as barcode_file, open('bowtie_input.fastq','w') as output_file:
	while True:
	    placeholder+=1
	    a1=allele_file.readline()
	    if a1 == '':
		break
            a2=allele_file.readline()
            a3=allele_file.readline()
            a4=allele_file.readline()
            b1=barcode_file.readline()
            b2=barcode_file.readline()
            b3=barcode_file.readline()
            b4=barcode_file.readline()
            quality_score=calc_quality(b4[:-1])
	    if quality_score >= 0.95: # if the barcode has <5% chance of error
		matches+=1 # write it to the bowtie input file
		output_file.write(a1[:-1]+'-'+b2)
		output_file.write(a2)
		output_file.write(a3)
		output_file.write(a4)

reverse_complement_file()
merge()
pair()
