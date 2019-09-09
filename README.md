# Barcode-Association

These scripts work together with Bowtie2 to associate barcodes with a-synuclein variants from long-read sequencing. They provide the variant dictionary used for the experiments described in [citation forthcoming]

Procedure
1) Run merge.py
2) Run Bowtie2 on the output of merge.py ('bowtie_input.fastq') against the list of reference sequences ('reference_sequences.fasta')
3) Run parse_bowtie.py

Associated Files
1) A list of well-represent barcodes ('barcode_list.txt')
2) A FASTA file of expected mutant sequences ('reference_sequences.fasta')
3) Raw sequencing data, available at NCBI [accession data forthcoming]
