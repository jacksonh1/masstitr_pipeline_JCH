'''
Run script with python2
I haven't added much documentation to this b/c it's a pretty short script
'''


import subprocess
import glob, os
import multiprocessing
import sys
import pandas as pd
import masstitr_tools as mt

# TODO use argparser instead of positional arguments


def run(file, reformat_command):
    # run BBduk reformat command to quality filter and de-interleave paired reads
    output_path = os.path.dirname(file)
    subprocess.call("{} in={} out1={}_f.fq out2={}_r.fq minavgquality=20".format(reformat_command, file, file, file), shell=True)
    subprocess.call("rm {}_r.fq".format(file), shell=True)
    # convert fastq to plain list of seqs (for Venkat's script)
    subprocess.call('python ./src/fq2str.py {}_f.fq'.format(file), shell=True)
    subprocess.call("rm {}_f.fq".format(file), shell=True)
    # run Venkat's script
    subprocess.call('python ./src/count_sequences.py {}_f_nts_only "{}/seq_counts" -c "{}/seq_counts_complete"'.format(file,output_path,output_path), shell=True)


def main():
    experiment_directory = str(sys.argv[1])
    sample_key_file = str(sys.argv[2])
    reformat_command = str(sys.argv[3])
    exp_dirs = mt.get_exp_dir_dict(sample_key_file, experiment_directory)
    for dir in exp_dirs.values():
        p = multiprocessing.Pool()
        for f in glob.glob(os.path.join(dir,"barcode*")):
            # skip files ending in ".fq". Makes it easier if you rerun this script
            if f.endswith(".fq"):
                continue
            if f.endswith("_nts_only"):
                continue
            # launch a process for each file (ish).
            # The result will be approximately one process per CPU core available.
            p.apply_async(run, [f, reformat_command])
        p.close()
        p.join() # Wait for all child processes to close.


if __name__ == "__main__":
    main()
