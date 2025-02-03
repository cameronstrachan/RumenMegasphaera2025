import os
from concurrent.futures import ThreadPoolExecutor
from functions import *

renamed_genomes_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/drep_select_mags_renamed/'

concatenated_genomes_file = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/drep_select_mags_renamed/concatenated_genomes.fa'

cpu_num = 36
blastn_cpu_num = 5

blast_db_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/blast_db_concatenated_genomes/'

if not os.path.exists(blast_db_dir):
    os.makedirs(blast_db_dir)

## run prodigal on concatenated genomes

def run_prodigal(input_file, output_dir):
        
    output_file_aa = f"{output_dir}{os.path.basename(input_file).split('.')[0]}.faa"
    output_file_nucl = f"{output_dir}{os.path.basename(input_file).split('.')[0]}.fna"
    prodigal_command = f"prodigal -p meta -i {input_file} -a {output_file_aa} -d {output_file_nucl}"
    
    if not os.path.exists(output_file_nucl):
        os.system(prodigal_command)
    
run_prodigal(concatenated_genomes_file, blast_db_dir)

## Run makeblastdb

prodigal_output_nucl = f"{blast_db_dir}{os.path.basename(concatenated_genomes_file).split('.')[0]}.fna"

makeblastdb_command = f"makeblastdb -in {prodigal_output_nucl} -dbtype nucl -out {blast_db_dir}concatenated_genomes"

if not os.path.exists(f"{blast_db_dir}concatenated_genomes.nhr"):
    os.system(makeblastdb_command)

## Blast mapped reads

mapped_reads_fasta_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/filtered_fasta/'

blast_output_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/blast_output_mapped_reads/'

if not os.path.exists(blast_output_dir):
    os.makedirs(blast_output_dir)

def blast_mapped_reads(fasta_file, output_dir, blast_db, blast_cpu_num=blastn_cpu_num):
            
    blast_output_file = f"{output_dir}{os.path.basename(fasta_file).split('.')[0]}.out"
    blast_command = f"blastn -db {blast_db} -query {fasta_file} -out {blast_output_file} -num_threads {blast_cpu_num} -max_target_seqs 1 -max_hsps 1 -outfmt '6 qseqid sseqid pident sstart send qstart qend evalue bitscore score qlen length'"
    
    if not os.path.exists(blast_output_file):
        os.system(blast_command)
        print(blast_command)

input_fasta_files = [f"{mapped_reads_fasta_dir}{file}" for file in os.listdir(mapped_reads_fasta_dir)]

with ThreadPoolExecutor(max_workers=cpu_num) as executor:
    executor.map(blast_mapped_reads, input_fasta_files, [blast_output_dir]*len(input_fasta_files), [f"{blast_db_dir}concatenated_genomes"]*len(input_fasta_files))
