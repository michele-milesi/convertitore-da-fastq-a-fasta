import Bio
import numpy as np
from Bio import Seq
from Bio.Seq import Seq
from Bio import SeqRecord
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

#Controlli su input?
L1 = int(input("L1: "))
L2 = int(input("L2: "))
Q1 = int(input("Q1: "))
Q2 = int(input("Q2: "))
P = float(input("P: "))
file_input_name = input("File to convert: ")
fastq_records = SeqIO.parse(file_input_name, 'fastq')
fastq_record_list = list(fastq_records)
out = open('./documents/michele/bioinformatica/assignment 2/output.fa', 'w')

for fastq_record in fastq_record_list:
    record_len = len(fastq_record.seq)
    phred_quality = fastq_record.letter_annotations['phred_quality']
    min_quality = min(phred_quality)
    intervals = []
    if record_len >= L1 and min_quality >= Q1 and record_len <= L2:
        bool_list = [quality >= Q2 for quality in phred_quality]
        start_list = [i for i in range(len(bool_list)) if bool_list[i] and (i == 0 or not bool_list[i - 1])]
        start_list[:0] = [1]
        end_list = [i for i in range(len(bool_list)) if bool_list[i] and (i == len(bool_list) - 1 or not bool_list[i + 1])]
        end_list[:0] = [0]
        intervals_length = [end_list[i] - start_list[i] + 1 for i in range(len(start_list))]
        max_breadth = max(intervals_length)
        index = intervals_length.index(max_breadth)
        start = start_list[index]
        end = end_list[index]
        if max_breadth >= P * record_len:   #   /100 (?)
            description = "len:" + str(record_len) + " min_quality:" + str(min_quality) + " "
            description += str(start_list[index]) + ":" + str(end_list[index])
            mean = sum(phred_quality[start:end + 1])/max_breadth if start <= end else float('NaN')
            description += " medium_quality:" + str(mean)
            fastq_record.description = description
            fasta_record = fastq_record.format('fasta')
            out.write(fasta_record)

out.close()
