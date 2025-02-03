for file in *.fastq.gz; do
    mv "$file" "${file//_1/_1_L001_R1_001}"
done

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path /data/Unit_LMM/selberherr-group/strachan/mega/amplicons/mcgloughlin2020/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /home/strachan/master/mega/amplicons/qiime_int_files/demux-single-end-mcgloughlin2020.qza

qiime demux summarize \
  --i-data /home/strachan/master/mega/amplicons/qiime_int_files/demux-single-end-mcgloughlin2020.qza \
  --o-visualization /home/strachan/master/mega/amplicons/qiime_int_files/demux-single-end-mcgloughlin2020.qzv
  
scp -r -i /Users/cameronstrachan/ssh/key_vetmed -P 12121 strachan@i121srv02.vu-wien.ac.at:/home/strachan/master/mega/amplicons/qiime_int_files/demux-single-end-mcgloughlin2020.qzv

qiime dada2 denoise-single \
  --i-demultiplexed-seqs /home/strachan/master/mega/amplicons/qiime_int_files/demux-single-end-mcgloughlin2020.qza \
  --p-trim-left 20 \
  --p-trunc-len 240 \
  --o-representative-sequences /home/strachan/master/mega/amplicons/qiime_int_files/rep-seqs-dada2-mcgloughlin2020.qza \
  --o-table /home/strachan/master/mega/amplicons/qiime_int_files/table-dada2-mcgloughlin2020.qza \
  --o-denoising-stats /home/strachan/master/mega/amplicons/qiime_int_files/stats-dada2-mcgloughlin2020.qza \
  --p-n-threads 60

qiime tools export \
  --input-path /home/strachan/master/mega/amplicons/qiime_int_files/table-dada2-mcgloughlin2020.qza \
  --output-path /home/strachan/master/mega/amplicons/qiime_int_files

biom convert -i /home/strachan/master/mega/amplicons/qiime_int_files/feature-table.biom -o /home/strachan/master/mega/amplicons/output/feature-table-100-mcgloughlin2020.txt --to-tsv

qiime tools export \
  --input-path /home/strachan/master/mega/amplicons/qiime_int_files/rep-seqs-dada2-mcgloughlin2020.qza \
  --output-path /home/strachan/master/mega/amplicons/output

mv /home/strachan/master/mega/amplicons/output/dna-sequences.fasta /home/strachan/master/mega/amplicons/output/dna-sequences-mcgloughlin2020.fasta
