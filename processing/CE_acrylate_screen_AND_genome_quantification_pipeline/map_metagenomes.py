import os
from concurrent.futures import ThreadPoolExecutor
from functions import *

ref_genomes_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/drep_select_mags/'
renamed_genomes_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/drep_select_mags_renamed/'

cpu_num = 180
input_file_ext = '.fa'

## Rename reference genomes

if not os.path.exists(renamed_genomes_dir):
    os.makedirs(renamed_genomes_dir)

def rename_ref_genomes(in_folder_path, rename_folder, threads=cpu_num, ext = input_file_ext):

    
    files = [file for file in os.listdir(in_folder_path) if file.endswith(ext)]
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        executor.map(rename_contigs, files, [in_folder_path]*len(files), [rename_folder]*len(files))


rename_ref_genomes(ref_genomes_dir, renamed_genomes_dir, cpu_num, input_file_ext)

## Concatenated reference genomes into one file

concatenated_genomes_file = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/drep_select_mags_renamed/concatenated_genomes.fa'

if not os.path.exists(concatenated_genomes_file):
    concatenate_command = f"cat {renamed_genomes_dir}*{input_file_ext} > {concatenated_genomes_file}"
    os.system(concatenate_command)

## Run BWA index
## Map metagenomes to reference genomes

metagenomes_1_dir = '/data/Unit_LMM/selberherr-group/strachan/malmuthuge_metagenomes/'
metagenomes_2_dir = '/data/Unit_LMM/selberherr-group/strachan/stewart_metagenomes_sub/'

bwa_index_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/bwa_index/'
bwa_output_dir = '/data/Unit_LMM/selberherr-group/strachan/rumen_mags/bwa_output/'

if not os.path.exists(bwa_index_dir):
    os.makedirs(bwa_index_dir)

if not os.path.exists(bwa_output_dir):
    os.makedirs(bwa_output_dir)

bwa_index_command = f"bwa index {concatenated_genomes_file} -p {bwa_index_dir}concatenated_genomes"

if not os.path.exists(f"{bwa_index_dir}concatenated_genomes.amb"):
    os.system(bwa_index_command)

mapping_database_path = f"{bwa_index_dir}concatenated_genomes"

list_fastq_file_path1 = [f"{metagenomes_1_dir}{file}" for file in os.listdir(metagenomes_1_dir) if file.endswith('.fastq')]
list_fastq_file_path2 = [f"{metagenomes_2_dir}{file}" for file in os.listdir(metagenomes_2_dir) if file.endswith('.fastq')]
list_fastq_file_path = list_fastq_file_path1 + list_fastq_file_path2

for fastq_file in list_fastq_file_path:
    
    fastq_file_name = os.path.basename(fastq_file).split('.')[0]
    sam_output_file = f"{bwa_output_dir}{fastq_file_name}.sam"
    bwa_map_command = f"bwa mem -t {cpu_num} {mapping_database_path} {fastq_file} > {sam_output_file}"
    
    if not os.path.exists(sam_output_file):
        os.system(bwa_map_command)

