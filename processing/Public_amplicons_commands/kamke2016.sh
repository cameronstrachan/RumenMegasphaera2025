for file in *.fastq.gz; do
    mv "$file" "${file//_1/_1_L001_R1_001}"
done

###

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/kamke2016/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/demux-single-end-kamke2016.qza

qiime demux summarize \
  --i-data /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/demux-single-end-kamke2016.qza \
  --o-visualization /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/demux-single-end-kamke2016.qzv


scp -r -i /Users/cameronstrachan/ssh/key_vetmed -P 12121 strachan@i121srv02.vu-wien.ac.at:/data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/demux-single-end-kamke2016.qzv /Users/cameronstrachan/Desktop/TEXT/Manucripts/

qiime dada2 denoise-single \
  --i-demultiplexed-seqs /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/demux-single-end-kamke2016.qza \
  --p-trim-left 20 \
  --p-trunc-len 300 \
  --o-representative-sequences /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/rep-seqs-dada2-kamke2016.qza \
  --o-table /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/table-dada2-kamke2016.qza \
  --o-denoising-stats /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/stats-dada2-kamke2016.qza \
  --p-n-threads 30

qiime tools export \
  --input-path /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/table-dada2-kamke2016.qza \
  --output-path /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files

biom convert -i /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/feature-table.biom -o /home/strachan/master/mega/amplicons/output/feature-table-100-kamke2016.txt --to-tsv

rm /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/feature-table.biom

qiime tools export \
  --input-path /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/qiime_int_files/rep-seqs-dada2-kamke2016.qza \
  --output-path /home/strachan/master/mega/amplicons/output

mv /home/strachan/master/mega/amplicons/output/dna-sequences.fasta /home/strachan/master/mega/amplicons/output/dna-sequences-kamke2016.fasta
