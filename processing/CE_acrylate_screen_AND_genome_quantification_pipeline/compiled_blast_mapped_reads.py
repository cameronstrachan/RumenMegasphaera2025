
import os
from concurrent.futures import ThreadPoolExecutor
import logging
import pandas as pd
import shutil

blast_trim_pident_cutoff = 95
blast_trim_aident_cutoff = 70
num_threads = 80

input_blast_data = "/data/Unit_LMM/selberherr-group/strachan/rumen_mags/blast_output_mapped_reads/"
output_blast_data = "/data/Unit_LMM/selberherr-group/strachan/rumen_mags/blast_out_reads_test/compiled_trimmed_blast_mapped_reads.csv"

def read_blastout(file_name, in_folder_path):
    # If file not empty
    if os.stat(os.path.join(in_folder_path, file_name)).st_size == 0:
        logging.warning('File %s is empty. Skipping...', file_name)
        return pd.DataFrame()
    else:
        read_path = os.path.join(in_folder_path, file_name)
        # Explicitly specifying the delimiter as tab
        df = pd.read_csv(read_path, header=None, delimiter='\t', encoding='utf-8')
        # Add the file name as a column
        df['file'] = file_name
        return df

def compile_trim_blast_results(in_folder_path, trimmed_output_file_path, num_threads=num_threads, blast_out_ext = '.out'):
    
    if not os.path.exists(trimmed_output_file_path):
    
        # Get all the files in the input folder
        blast_out_files = [file for file in os.listdir(in_folder_path) if file.endswith(blast_out_ext)]

        # Concatenate all the files into a pandas dataframe
        logging.info('Concatenating files...')
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            dfs = list(executor.map(read_blastout, blast_out_files, [in_folder_path]*len(blast_out_files)))

        # Concatenate all DataFrames
        df = pd.concat(dfs, ignore_index=True)

        # Add column names (qseqid sseqid pident sstart send qstart qend evalue bitscore score qlen length)
        df.columns = ['qseqid', 'sseqid', 'pident', 'sstart', 'send', 'qstart', 'qend', 'evalue', 'bitscore', 'score', 'qlen', 'length', 'file']

        # Calculate alignment percentage
        df['aident'] = (df['length'] / df['qlen']) * 100

        # Select only columns we want (qseqid sseqid pident aident)
        df = df[['qseqid', 'sseqid', 'pident', 'aident', 'file']]

        logging.info('Trimming dataframe...')
        df = df[(df['pident'] > blast_trim_pident_cutoff) & (df['aident'] > blast_trim_aident_cutoff)]

        # Save the dataframe to a csv
        logging.info('Saving dataframe...')
        df.to_csv(trimmed_output_file_path, index=False)

    else:
        logging.info(f"File {trimmed_output_file_path} already exists, skipping.")

compile_trim_blast_results(input_blast_data, output_blast_data)

