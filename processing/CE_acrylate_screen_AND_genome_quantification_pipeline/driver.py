import os
from concurrent.futures import ThreadPoolExecutor
import logging
import pandas as pd
from functions import *
import shutil

####################################################################################################
# Driver script to rename genome files, run prodigal, and blastp a file against the orfs.
# Currently only works at aa level
# In the control Megasphaera genomes, the contig number is getting flipped, so gotta figure that out. 
####################################################################################################

# Setup basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

input_folder_path = 'mags/'
rename_folder_path = 'mags_renamed/'
prodigal_folder_path = 'mags_prodigal/'
db_folder_path = 'blast_db_aa/'
query_file_path = 'marker_proteins.fasta'
output_folder_path = 'blast_output/'

# create tmp dir if doesn't exist
if not os.path.exists('tmp'):
    os.makedirs('tmp')

compiled_prodigal_orf_positions_path = 'tmp/' + 'compiled_prodigal_orf_positions.csv'
trimmed_output_file_path = 'tmp/' + 'trimmed_blast_results.csv'


final_merged_file_path = 'merged_trimmed_blast_results_orf_positions.csv'

path_lst = [input_folder_path, rename_folder_path, prodigal_folder_path, db_folder_path, output_folder_path]

num_threads = 120
genome_ext = '.fa'
prodigal_ext = '.faa'
blast_out_ext = '.out'
blast_db_type = 'prot'

blast_trim_pident_cutoff = 40
blast_trim_aident_cutoff = 50

def predict_orfs_genomes(in_folder_path, rename_folder, out_folder_path, threads=4, ext = '.fa'):
    """
    Run renaming and prodigal on all files in a folder in parallel.
    """
    
    files = [file for file in os.listdir(in_folder_path) if file.endswith(ext)]
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        executor.map(rename_contigs, files, [in_folder_path]*len(files), [rename_folder]*len(files))

    with ThreadPoolExecutor(max_workers=threads) as executor:
        executor.map(run_prodigal, files, [rename_folder]*len(files), [out_folder_path]*len(files))

def compile_prodigal_orf_positions(in_folder_path, out_file_path, ext='.faa', num_threads=4):

    if not os.path.exists(out_file_path):

        # Get all the file paths
        prodigal_file_paths = [os.path.join(in_folder_path, file) for file in os.listdir(in_folder_path) if file.endswith(prodigal_ext)]

        # Use ThreadPoolExecutor to process files in parallel
        data = []
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            results = executor.map(parse_prodigal_file, prodigal_file_paths)
            for result in results:
                data.extend(result)

        # Create a DataFrame
        df = pd.DataFrame(data, columns=['sseqid', 'start', 'stop', 'direction'])

        # Save the DataFrame to a CSV file
        df.to_csv(out_file_path, index=False)

    else:
        logging.info(f"File {out_file_path} already exists, skipping.")

def blast_orfs(query_file_path, in_folder_path, db_folder_path, out_folder_path, db_type='prot', threads=num_threads, ext='.faa'):
    """
    Create blast databases from FASTA files and run blastp on a given file in parallel.
    """

    prodigal_files = [file for file in os.listdir(in_folder_path) if file.endswith(ext)]
    db_files = [file.split(ext)[0] for file in prodigal_files]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        executor.map(run_makeblastdb, prodigal_files, [in_folder_path]*len(prodigal_files), [db_folder_path]*len(prodigal_files), [db_type]*len(prodigal_files))

    with ThreadPoolExecutor(max_workers=threads) as executor:
        executor.map(run_blastp, [query_file_path]*len(db_files), [out_folder_path]*len(db_files), [db_folder_path]*len(db_files), db_files)

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
        df.columns = ['qseqid', 'sseqid', 'pident', 'sstart', 'send', 'qstart', 'qend', 'evalue', 'bitscore', 'score', 'qlen', 'length']

        # Calculate alignment percentage
        df['aident'] = (df['length'] / df['qlen']) * 100

        # Select only columns we want (qseqid sseqid pident aident)
        df = df[['qseqid', 'sseqid', 'pident', 'aident']]

        # Splitting the 'qseqid' column
        split_qseqid_df = df['sseqid'].str.rsplit('_', n=2, expand=True)

        # Naming the new columns
        split_qseqid_df.columns = ['genome', 'contig', 'orf']

        # Joining the new columns with the original DataFrame
        df = df.join(split_qseqid_df)

        # Trim the dataframe for pident > 50 and aident > 50
        logging.info('Trimming dataframe...')
        df = df[(df['pident'] > blast_trim_pident_cutoff) & (df['aident'] > blast_trim_aident_cutoff)]

        # Save the dataframe to a csv
        logging.info('Saving dataframe...')
        df.to_csv(trimmed_output_file_path, index=False)

    else:
        logging.info(f"File {trimmed_output_file_path} already exists, skipping.")

def merge_trimmed_blast_results_and_orf_positions(compiled_orf_positions_file_path, trimmed_blast_results_file_path, out_file_path):
    
    if not os.path.exists(out_file_path):

        # Read the CSV files
        df_orf_positions = pd.read_csv(compiled_orf_positions_file_path)
        df_blast_results = pd.read_csv(trimmed_blast_results_file_path)

        # Merge the DataFrames
        merged_df = pd.merge(df_blast_results, df_orf_positions, on='sseqid', how='left')

        # Save the merged DataFrame to a new CSV file
        merged_df.to_csv(out_file_path, index=False)

        print("Merging complete. Data saved to", out_file_path)

    else:
        logging.info(f"File {out_file_path} already exists, skipping.")


def main():
    """
    Main function to orchestrate the pipeline.
    """
    if not all(os.path.exists(path) for path in path_lst):
        logging.error("One or more specified directories do not exist.")
        return

    predict_orfs_genomes(input_folder_path, rename_folder_path, prodigal_folder_path, threads=num_threads, ext=genome_ext)

    blast_orfs(query_file_path, prodigal_folder_path, db_folder_path, output_folder_path, db_type=blast_db_type, threads=num_threads, ext=prodigal_ext)

    compile_prodigal_orf_positions(prodigal_folder_path, compiled_prodigal_orf_positions_path, num_threads=num_threads, ext=prodigal_ext)

    compile_trim_blast_results(output_folder_path, trimmed_output_file_path, num_threads=num_threads, blast_out_ext=blast_out_ext)

    merge_trimmed_blast_results_and_orf_positions(compiled_prodigal_orf_positions_path, trimmed_output_file_path, final_merged_file_path)

    # Remove tmp dir
    shutil.rmtree('tmp')

if __name__ == '__main__':
    main()
