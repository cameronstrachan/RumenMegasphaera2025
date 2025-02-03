import os
import subprocess
import logging
import pandas as pd

def rename_contigs(file_name, in_folder_path, out_folder_path):
    """
    Rename contigs in a FASTA file.
    """
    try:
        input_file_path = os.path.join(in_folder_path, file_name)
        output_file_path = os.path.join(out_folder_path, file_name)

        if not os.path.exists(output_file_path):

            with open(input_file_path, 'r') as f_in, open(output_file_path, 'w') as f_out:
                file_prefix = os.path.splitext(file_name)[0]
                contig_count = 1

                for line in f_in:
                    if line.startswith('>'):
                        new_header = f'>{file_prefix}_{contig_count}\n'
                        f_out.write(new_header)
                        contig_count += 1
                    else:
                        f_out.write(line)
        else:
            logging.info(f"File {output_file_path} already exists, skipping.")
    except Exception as e:
        logging.error(f"Error renaming file {file_name}: {e}")


def run_prodigal(file_name, in_folder_path, out_folder_path, prodigal_ext = '.faa'):
    """
    Run Prodigal on a given file.
    """
    try:
        file_prefix = os.path.splitext(file_name)[0]
        input_file_path = os.path.join(in_folder_path, file_name)
        output_file_path = os.path.join(out_folder_path, file_prefix + prodigal_ext)

        if not os.path.exists(output_file_path):
            command = f"prodigal -i {input_file_path} -a {output_file_path} -p meta"
            subprocess.run(command, shell=True)
        else:
            logging.info(f"File {output_file_path} already exists, skipping.")
    except Exception as e:
        logging.error(f"Error running Prodigal on file {file_name}: {e}")

def parse_prodigal_file(file_path):
    temp_data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                parts = line.split('#')
                sseqid = parts[0].split()[0][1:]  # Remove '>' and get the first part
                start = parts[1].strip()
                stop = parts[2].strip()
                direction = parts[3].strip()
                temp_data.append([sseqid, start, stop, direction])
    return temp_data


def run_makeblastdb(file_name, in_folder_path, out_folder_path, db_type='prot'):
    """
    Run makeblastdb on a given file.
    """
    try:
        file_prefix = os.path.splitext(file_name)[0]
        input_file_path = os.path.join(in_folder_path, file_name)
        output_file_path = os.path.join(out_folder_path, file_prefix)

        if db_type == 'prot':
            blastdb_ext = '.phr'
        else:
            blastdb_ext = '.nhr'

        if not os.path.exists(output_file_path + blastdb_ext):
            command = f"makeblastdb -in {input_file_path} -dbtype {db_type} -out {output_file_path}"
            subprocess.run(command, shell=True)
        else:
            logging.info(f"File {output_file_path + blastdb_ext} already exists, skipping.")
    except Exception as e:
        logging.error(f"Error running makeblastdb on file {file_name}: {e}")

def run_blastp(file_name, out_folder_path, db_folder, db_file):
    """
    Run blastp on a given file.
    """
    try:
        file_prefix = file_name.split(".f")[0]
        output_file_path = os.path.join(out_folder_path, f"{file_prefix}_vs_{db_file}.out")
        db_path = os.path.join(db_folder, db_file)

        if not os.path.exists(output_file_path):

            command = f"blastp -query {file_name} -db {db_path} -out {output_file_path} -outfmt '6 qseqid sseqid pident sstart send qstart qend evalue bitscore score qlen length' -evalue 1e-5"
            subprocess.run(command, shell=True)
        else:
            logging.info(f"File {output_file_path} already exists, skipping.")

    except Exception as e:
        logging.error(f"Error running blastp on file {file_name}: {e}")

def read_blastout(file_name, in_folder_path):
    # If file not empty
    if os.stat(os.path.join(in_folder_path, file_name)).st_size == 0:
        logging.warning('File %s is empty. Skipping...', file_name)
        return pd.DataFrame()
    else:
        read_path = os.path.join(in_folder_path, file_name)
        # Explicitly specifying the delimiter as tab
        return pd.read_csv(read_path, header=None, delimiter='\t', encoding='utf-8')