import os, sys
import re

# map metagenomes to selected genomes competitively

threads = 80
genomes_folder = '/data/Unit_LMM/selberherr-group/strachan/mega/transcriptomes/genomes/'
metagenome_reads_folder = '/data/Unit_LMM/selberherr-group/strachan/mega/transcriptomes/park2022/'
intermediate_folder = '/data/Unit_LMM/selberherr-group/strachan/mega/transcriptomes/bwa_intermediate_files/'
tables_folder = '/home/strachan/master/mega/transcriptomes/output_park_veillonella/'

ref = genomes_folder +  'veillonellales_mags.fasta'
ref_check = genomes_folder +  'bwa_index/veillonellales_mags.fasta.bwt'

if not os.path.exists(ref_check):
        command0 = 'bwa index ' + ref
        print(command0 + '\n')
        os.system(command0)
        os.system('mv ' + genomes_folder + 'veillonellales_mags.fasta.* ' + genomes_folder + 'bwa_index/')

fastq_files = [f for f in os.listdir(metagenome_reads_folder) if f.endswith(('_1.fastq', '_2.fastq'))]
sample_names = list(set([re.sub('_\d.fastq', '', s) for s in fastq_files]))

database = genomes_folder  + 'bwa_index/veillonellales_mags.fasta'

for sample in sample_names:
        
        forward_read_trim = metagenome_reads_folder + sample + '_1.fastq'
    
        output_file_loc = intermediate_folder + sample + '.R1.metagenome.veillonellales_mags.sam'
        output_file_loc_sam_sort = intermediate_folder + sample + '.R1.metagenome.veillonellales_mags.sort.sam'

        if not os.path.exists(output_file_loc_sam_sort):
                command1 = 'bwa mem -t ' + str(threads) + ' ' + database + ' ' + forward_read_trim  + ' > ' + output_file_loc
                print(command1 + '\n')
                os.system(command1)
                
                command2 = 'samtools sort -O sam -@ ' + str(threads) + ' ' + output_file_loc + ' -o ' + output_file_loc_sam_sort
                print(command2 + '\n')
                os.system(command2)

prokka_folder = '/data/Unit_LMM/selberherr-group/strachan/mega/transcriptomes/prokka_ref/'


sam_files = [f for f in os.listdir(intermediate_folder) if f.endswith('.R1.metagenome.veillonellales_mags.sort.sam')]

for sam_file in sam_files:
        output_table = tables_folder + sam_file.split('.sa')[0] + '.txt'
        
        if not os.path.exists(output_table):
                command3 = 'htseq-count -s no -t CDS -i ID --additional-attr=gene --additional-attr=product ' + intermediate_folder + sam_file + ' ' + prokka_folder + 'veillonellales_mags.gff' + ' > ' + output_table 
                print(command3 + '\n')
                os.system(command3)
