for file in *.fastq.gz; do
    mv "$file" "${file//_1/_1_L001_R1_001}"
done

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path /data/Unit_LMM/selberherr-group/strachan/amplicons/gaire2023/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/demux-single-end-gaire2023.qza

qiime demux summarize \
  --i-data /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/demux-single-end-gaire2023.qza \
  --o-visualization /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/demux-single-end-gaire2023.qzv
  
scp -r -i /Users/cameronstrachan/ssh/key_vetmed -P 12121 strachan@i121srv02.vu-wien.ac.at:/data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/demux-single-end-gaire2023.qzv 

qiime dada2 denoise-single \
  --i-demultiplexed-seqs /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/demux-single-end-gaire2023.qza \
  --p-trim-left 20 \
  --p-trunc-len 240 \
  --o-representative-sequences /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/rep-seqs-dada2-gaire2023.qza \
  --o-table /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/table-dada2-gaire2023.qza \
  --o-denoising-stats /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/stats-dada2-gaire2023.qza \
  --p-n-threads 50

qiime tools export \
  --input-path /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/table-dada2-gaire2023.qza \
  --output-path /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files

biom convert -i /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/feature-table.biom -o /home/strachan/master/amplicons/output/feature-table-100-gaire2023.txt --to-tsv

qiime tools export \
  --input-path /data/Unit_LMM/selberherr-group/strachan/amplicons/qiime_int_files/rep-seqs-dada2-gaire2023.qza \
  --output-path /home/strachan/master/amplicons/output

mv /home/strachan/master/amplicons/output/dna-sequences.fasta /home/strachan/master/amplicons/output/dna-sequences-gaire2023.fasta
