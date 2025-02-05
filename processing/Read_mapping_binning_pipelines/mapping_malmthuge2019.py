import os, sys
import re

# map metagenomes to selected genomes competitively

threads = 120
threads_sam = 20

genomes_folder = '/data/Unit_LMM/selberherr-group/strachan/select4_genomes/'
metagenome_reads_folder = '/data/Unit_LMM/selberherr-group/strachan/malmuthuge_metagenomes/'
intermediate_folder = '/data/Unit_LMM/selberherr-group/strachan/select4_bwa_intermediate_files/'

ref = genomes_folder +  'Megasphaera_genomes.fasta'
ref_check = genomes_folder +  'bwa_index/Megasphaera_genomes.fasta.amb'

if not os.path.exists(ref_check):
        command0 = 'bwa index ' + ref
        print(command0 + '\n')
        os.system(command0)
        os.system('mv ' + genomes_folder + 'Megasphaera_genomes.fasta.* ' + genomes_folder + 'bwa_index/')

fastq_files = [f for f in os.listdir(metagenome_reads_folder) if f.endswith('_1.fastq')]
sample_names = list(set([re.sub('_\d.fastq', '', s) for s in fastq_files]))

database = genomes_folder  + 'bwa_index/Megasphaera_genomes.fasta'

for sample in sample_names:
        
        forward_read_trim = metagenome_reads_folder + sample + '_1.fastq'

        output_file_loc = intermediate_folder + sample + '.R1.metagenome.select4_ref_genomes.sam'
        output_file_loc_sam_mapped = intermediate_folder + sample + '.R1.metagenome.select4_ref_genomes.mapped.bam'
        output_file_loc_sam_mapped_sort = intermediate_folder + sample + '.R1.metagenome.select4_ref_genomes.mapped.sort.sam'

        if not os.path.exists(output_file_loc):
                command1 = 'bwa mem -t ' + str(threads) + ' ' + database + ' ' + forward_read_trim  + ' > ' + output_file_loc
                print(command1 + '\n')
                os.system(command1)
                
                command2 = 'samtools view -b -F 4 -@ ' + str(threads_sam) + ' ' + output_file_loc + ' -o ' + output_file_loc_sam_mapped
                print(command2 + '\n')
                os.system(command2)

                command3 = 'samtools sort -O sam -@ ' + str(threads_sam) + ' ' + output_file_loc_sam_mapped + ' -o ' + output_file_loc_sam_mapped_sort
                print(command3 + '\n')
                os.system(command3)

                os.system('rm ' + output_file_loc)
