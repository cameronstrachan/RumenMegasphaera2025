import os
from concurrent.futures import ThreadPoolExecutor
from functions import *

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

cpu_num = 180
input_file_ext = '.fa'

bwa_output_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/bwa_output/'

## Filter SAM for mapped reads

filtered_sam_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/filtered_sam/'

if not os.path.exists(filtered_sam_dir):
    os.makedirs(filtered_sam_dir)

def filter_sam(sam_file, output_dir):
        
        filtered_sam_file = f"{output_dir}{os.path.basename(sam_file)}"
        filter_command = f"samtools view -F 4 {sam_file} > {filtered_sam_file}"
        
        if not os.path.exists(filtered_sam_file):
            os.system(filter_command)

input_sam_files = [f"{bwa_output_dir}{file}" for file in os.listdir(bwa_output_dir) if file.endswith('.sam')]

with ThreadPoolExecutor(max_workers=cpu_num) as executor:
    executor.map(filter_sam, input_sam_files, [filtered_sam_dir]*len(input_sam_files))

## Convert SAM to fasta

filtered_fasta_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/filtered_fasta/'

if not os.path.exists(filtered_fasta_dir):
    os.makedirs(filtered_fasta_dir)

def sam_to_fasta_biopython(sam_file, output_dir):
    
    fasta_file = f"{output_dir}{os.path.basename(sam_file).replace('.sam', '.fasta')}"
    
    with open(sam_file) as samfile:
        with open(fasta_file, 'w') as fastafile:
            for line in samfile:
                if not line.startswith('@'):
                    line = line.split('\t')
                    read_id = line[0]
                    read_seq = line[9]
                    fastafile.write(f">{read_id}\n{read_seq}\n")

input_sam_files = [f"{filtered_sam_dir}{file}" for file in os.listdir(filtered_sam_dir) if file.endswith('.sam')]

with ThreadPoolExecutor(max_workers=cpu_num) as executor:
    executor.map(sam_to_fasta_biopython, input_sam_files, [filtered_fasta_dir]*len(input_sam_files))
