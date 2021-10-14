import pandas as pd
import subprocess
import sys
import glob, os

# TODO use argparser instead of positional arguments


def organize_files(barcode_directory, experiment_directory, exp, barcodes):
    output_dir = os.path.join(experiment_directory, exp)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    for bcode in barcodes:
        in_file = os.path.join(barcode_directory, bcode)
        out_file = os.path.join(output_dir, bcode)
        subprocess.call("mv {} {}".format(in_file, out_file), shell=True)


def main():
    barcode_directory = str(sys.argv[1])
    experiment_directory = str(sys.argv[2])
    sample_key_file = str(sys.argv[3])
    samples = pd.read_csv(sample_key_file)
    experiments = list(samples['experiment'].unique())
    for exp in experiments:
        barcodes = samples[samples['experiment']==exp]['barcode'].to_list()
        organize_files(barcode_directory, experiment_directory, exp, barcodes)


if __name__ == "__main__":
    main()
