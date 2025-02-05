# python libraries
import os, sys
import subprocess

# custom libraries
system = str(input('\n' + 'Local or Server (L or S):'))

if system == 'S':
    sys.path.insert(0, '/home/strachan/master/')
else:
    sys.path.insert(0, '/Users/cameronstrachan/master/')

#bwa index dataflow/04-coassembly/final.contigs.fa

files = [f for f in os.listdir('dataflow/01-fastq/') if f.endswith("_R1.fastq.gz")]

for file in files:

    file_prefix = file.split('_R1')[0]

    file_r1 = file_prefix + '_trimmed_R1.fastq.gz'
    file_r2 = file_prefix + '_trimmed_R2.fastq.gz'

    sam_file = file_prefix + '_coassembly_' + ".sam"
    command = "bwa mem -t 70 dataflow/04-coassembly/final.contigs.fa " + 'dataflow/01-fastq/' + file_r1  + ' dataflow/01-fastq/' + file_r2  + " > " + "dataflow/03-alignments/" + sam_file

    os.system(command)

    bam_file = file_prefix + '_coassembly_' + ".bam"
    command = "samtools view -@ 70 -bS " + "dataflow/03-alignments/" + sam_file + " > " + "dataflow/03-alignments/" + bam_file
    os.system(command)

    bam_file_sorted = file_prefix + '_coassembly_' + ".sorted.bam"
    command = "samtools sort -@ 70 " + "dataflow/03-alignments/" + sam_file + " > " + "dataflow/03-alignments/" + bam_file_sorted
    os.system(command)

command = 'jgi_summarize_bam_contig_depths --outputDepth dataflow/04-tables/megahit-metagenome-depth.txt dataflow/03-alignments/*sorted.bam'
os.system(command)

command = 'metabat2 -i dataflow/04-coassembly/final.contigs.fa -a dataflow/04-tables/megahit-metagenome-depth.txt -o dataflow/04-bins/bin -m 1500 -t 70 -v --unbinned'
os.system(command)
