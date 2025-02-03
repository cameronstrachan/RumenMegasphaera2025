import os
import csv

# Hardcoded directories
directory1 = '/data/Unit_LMM/selberherr-group/strachan/malmuthuge_metagenomes'
directory2 = '/data/Unit_LMM/selberherr-group/strachan/stewart_metagenomes_sub'

# Output CSV file name
output_csv = 'total_read_counts.csv'

def count_reads_in_fastq(file_path):
    """
    Counts the number of reads in a FASTQ file.
    Each read in a FASTQ file consists of 4 lines.
    """
    try:
        with open(file_path, 'r') as file:
            line_count = sum(1 for _ in file)
        read_count = line_count // 4
        return read_count
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return 0

def find_fastq_files(directory):
    """
    Recursively finds all FASTQ files in the given directory.
    """
    fastq_files = []
    for root, _, files in os.walk(directory):
        for name in files:
            if name.endswith(('.fastq', '.fq')):
                fastq_files.append(os.path.join(root, name))
    return fastq_files

# Main script
if __name__ == '__main__':
    # Find all FASTQ files in both directories
    fastq_files = find_fastq_files(directory1) + find_fastq_files(directory2)

    # Count reads in each file and store the results
    results = []
    for file_path in fastq_files:
        read_count = count_reads_in_fastq(file_path)
        results.append({'Filename': os.path.basename(file_path), 'ReadCount': read_count})

    # Write results to a CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Filename', 'ReadCount']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for entry in results:
            writer.writerow(entry)

    print(f"Read counts have been written to {output_csv}")
