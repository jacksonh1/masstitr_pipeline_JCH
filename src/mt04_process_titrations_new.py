import numpy as np
import pandas as pd
import glob
import os
import numpy as np
import sys
import masstitr_tools as mt
from Bio import Seq

# %%
# ==============================================================
# // parameters
# ==============================================================

# count cutoff threshold
FILTER_METHOD='cells'
COUNT_CUTOFF = 10
CLUSTER_READS = True
# list of barcodes to exclude from the analysis at this stage
EXCLUDED_BARCODES = ["barcode_5", "barcode_20", "barcode_21"]

# weighting of the different bins
WEIGHTS = {
    1:3000,
    2:300,
    3:30,
    4:0
}



# %%
# ==============================================================
# // group gates by concentrations, calculate binding signal for each concentration
# ==============================================================


def calc_binding_signal(df, conc, bin_key, weights):
    df[str(conc)]=(df[bin_key[1]]*weights[1] + df[bin_key[2]]*weights[2] + df[bin_key[3]]*weights[3] + df[bin_key[4]]*weights[4])/(df[bin_key[1]] + df[bin_key[2]] + df[bin_key[3]] + df[bin_key[4]])
    return df[['seq'] + [str(conc)]]


def process_exp(exp_directory, experiment, sample_key, count_cutoff, excluded_barcodes, weights):
    '''run pipeline for a single experiment'''
    file_list = [f for f in glob.glob(os.path.join(exp_directory, "seq_counts_complete", "*nts_only"))]
    R, T, barcode_cols = mt.load_and_merge_data(file_list, excluded_barcodes)
    R['AA_seq'] = R['seq'].apply(mt.trans_string)
    if CLUSTER_READS:
        R = mt.cluster_AA_seqs(R, barcode_cols)
    read_f = mt.read_count2fraction(R, T, barcode_cols)
    # get cell counts
    x = mt.read_fraction2cell_count(read_f, sample_key, barcode_cols)
    # run main processing
    b_curve = processed_data2binding_curve(x, R, weights, sample_key, experiment, count_cutoff, n_missing_allowed=2)
    return R, x, b_curve


def processed_data2binding_curve(x, R, weights, sample_key, experiment, count_cutoff, n_missing_allowed=2):
    '''calculate binding signal for a single experiment - used in `process_exp()`'''
    concentrations = sorted(list(sample_key['concentration (nM)'].unique()), reverse=True)
    b_curve = pd.DataFrame(index=x['seq'])
    for conc in concentrations:

        # get list of barcodes for gates used at `conc`
        barcodes = mt.conc2barcodes(sample_key, experiment, conc)

        # get bin # for each barcode. Bin # needs to correspond to
        # numbering in `weights`
        bin_key = mt.barcodes2bin(sample_key, barcodes)

        # select data for only the gates used at `conc`
        x_conc = x[['seq'] + barcodes].copy()
        R_conc = R[['seq'] + barcodes].copy()

        # filter
        if FILTER_METHOD=='reads':
            x_conc = mt.single_conc_readcount_filter(x_conc, R_conc, count_cutoff)
        if FILTER_METHOD=='cells':
            x_conc = mt.single_conc_cell_count_filter(x_conc, count_cutoff)

        #calculate binding signal
        df = calc_binding_signal(
            x_conc,
            conc,
            bin_key,
            weights
        )

        # merge to create full curve
        # 'outer' merge will keep all sequences even if for some
        # concentrations they have <count_cutoff
        b_curve = b_curve.merge(df, on='seq', how='outer')
    b_curve = b_curve[b_curve.isna().sum(1)<=n_missing_allowed]
    return b_curve


def save_exp_params(output_file):
    with open(output_file, 'w') as handle:
        handle.write(f"{FILTER_METHOD=}\n")
        handle.write(f"{COUNT_CUTOFF=}\n")
        handle.write(f"{CLUSTER_READS=}\n")
        handle.write(f"{EXCLUDED_BARCODES=}\n")
        handle.write(f"{WEIGHTS=}\n")


def main(count_cutoff, excluded_barcodes, weights):
    experiment_directory = str(sys.argv[1])
    sample_key_file = str(sys.argv[2])
    binding_curve_folder = str(sys.argv[3])

    sample_key = pd.read_csv(sample_key_file)
    exp_dirs = mt.get_exp_dir_dict(sample_key_file, experiment_directory)
    if not os.path.exists(binding_curve_folder):
        os.mkdir(binding_curve_folder)
    for exp, exp_directory in exp_dirs.items():
        print("experiment: {}".format(exp_directory))
        R, x, b_curve = process_exp(
            exp_directory=exp_directory,
            experiment=exp,
            sample_key=sample_key,
            count_cutoff=count_cutoff,
            excluded_barcodes=excluded_barcodes,
            weights=weights
        )
        # x.to_csv('{}{}-cell_counts.csv'.format(
        #         os.path.join(binding_curve_folder,''),
        #         exp
        #     ), index=False)
        b_curve.to_csv('{}{}-binding_curves.csv'.format(
                os.path.join(binding_curve_folder,''),
                exp
            ), index=False)
        # R.to_csv('{}{}-read_counts.csv'.format(
        #         os.path.join(binding_curve_folder,''),
        #         exp
        #     ), index=False)
        save_exp_params(os.path.join(binding_curve_folder,'mt04_processing_parameters.txt'))


if __name__ == "__main__":
    main(COUNT_CUTOFF, EXCLUDED_BARCODES, WEIGHTS)


