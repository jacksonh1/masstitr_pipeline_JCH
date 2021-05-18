import numpy as np
import pandas as pd
import glob
import os
import numpy as np
import sys

# %%
# ==============================================================
# // parameters
# ==============================================================

# read count cutoff threshold
# count_cutoff = 100
count_cutoff = 1

# list of barcodes to exclude from the analysis at this stage
excluded_barcodes = ["barcode_5", "barcode_20", "barcode_21"]

# weighting of the different bins
weights = {
    1:3000,
    2:300,
    3:30,
    4:0
}

# %%
# ==============================================================
# // functions
# ==============================================================


# import data --------------------------------------------------
# TODO: add documentation
def df_import_1(file):
    base = os.path.basename(file)
    basenoext = os.path.splitext(base)[0].replace('_f_nts_only', '')
    df = pd.read_csv(
            file,
            sep='\t',
            header=None,
            names=[basenoext, 'useless', 'seq']
        ).drop("useless", axis=1)
    return df


def load_and_merge_data(file_list, excluded_barcodes):
    """
    TODO: output from: {script name}
    load NGS data (output from: ) and merge into a single DataFrame
    """
    R = pd.DataFrame(columns=['seq'])
    for f in file_list:
        if os.path.splitext(os.path.basename(f))[0].replace('_f_nts_only', '') in excluded_barcodes:
            continue
        df1 = df_import_1(f)
        R = pd.merge(R, df1, on='seq', how='outer')
    R = R.fillna(0)
    # reorder columns
    barcode_cols = [col for col in R.columns]
    barcode_cols.sort()
    barcode_cols.remove('seq')
    R = R[['seq'] + barcode_cols]
    T = R[barcode_cols].sum()
    return R, T, barcode_cols


def initial_filter(R, count_cutoff):
    '''
    filter out any sequences where all of the bins are below cutoff threshold
    i.e. keep any sequences where at least one of the barcodes has >100 reads
    '''
    return R[(R.drop("seq", axis=1) > count_cutoff).any(1)]


# normalization --------------------------------------------------
def read_count2fraction(R_filtered, T, barcode_cols):
    '''
    normalize by dividing read counts for each sequence by total number of
    reads in the sample.
    '''
    # read fraction: R/T
    read_f = R_filtered.copy() # probably unnecessary but let's keep the original dataframe anyway
    for col in barcode_cols:
        read_f[col] = R_filtered[col]/T[col]
    return read_f


def read_fraction2cell_count(read_f, sample_key, barcode_cols):
    """
    For each sample, multiply the read fractions by the number of cells in the bin
        x = cells*(R/T)
    """
    # multiply by cells: x = cells*(R/T)
    x = read_f.copy()
    for col in barcode_cols:
        c_ijk = sample_key[sample_key['barcode']==col]['cell_count'].values[0]
        x[col] = x[col]*c_ijk
    return x


# constucting binding curve---------------------------------------
def conc2barcodes(sample_key, exp, conc):
    return list(sample_key[(sample_key['experiment']==exp) & (sample_key['concentration (nM)']==conc)]['barcode'])


def barcodes2bin(sample_key, barcodes):
    bin_key = {}
    for barcode in barcodes:
        b = sample_key[sample_key['barcode']==barcode]['bin'].values[0]
        bin_key[b] = barcode
    return bin_key


def single_conc_readcount_filter(x_df, R_df, count_cutoff):
    """filter out sequences that don't have at least `count_cutoff` reads in one of the barcodes of `x_df`

    Parameters
    ----------
    x_df : pandas DataFrame
        DataFrame of the normalized/processed data (x) corresponding to
        the gates in a single concentration point. `x_df` should also
        have a column of the sequences labelled "seq".
        This is usually a subset of the full x DataFrame
        example: `x[['seqs', 'barcode_168', 'barcode_169', ... ]]
    R_df : pandas DataFrame
        DataFrame of the raw read counts corresponding to
        the gates in a single concentration point. `R_df` should also
        have a column of the sequences labelled "seq". The columns
        should correspond to the columns in `x_df`
        This is usually a subset of the full R DataFrame
        example: `R[['seqs', 'barcode_168', 'barcode_169', ... ]]
    count_cutoff : int
        The read count cutoff to filter the sequences by

    Returns
    -------
    pandas DataFrame
        dataframe containing only the sequences passing the filter
    """

    passing_sequences = R_df[(R_df.drop("seq", axis=1) >= count_cutoff).any(1)].seq.to_list()
    x_df = x_df[x_df['seq'].isin(passing_sequences)]
    return x_df


def calc_binding_signal(df, conc, barcodes, bin_key, weights):
    df[str(conc)]=(df[bin_key[1]]*weights[1] + df[bin_key[2]]*weights[2] + df[bin_key[3]]*weights[3] + df[bin_key[4]]*weights[4])/(df[bin_key[1]] + df[bin_key[2]] + df[bin_key[3]] + df[bin_key[4]])
    return df[['seq'] + [str(conc)]]


def processed_data2binding_curve(x, R_filtered, weights, sample_key, experiment, count_cutoff, n_missing_allowed=2):
    concentrations = list(sample_key['concentration (nM)'].unique())
    b_curve = pd.DataFrame(index=x['seq'])
    for conc in concentrations:
        # get list of barcodes for concentration
        barcodes = conc2barcodes(sample_key, experiment, conc)
        # get bin # for each barcode. Bin # needs to correspond to
        # numbering in `weights`
        bin_key = barcodes2bin(sample_key, barcodes)
        # select data for only the gates used at `conc`
        x_conc = x[['seq'] + barcodes].copy()
        R_filt_conc = R_filtered[['seq'] + barcodes].copy()
        # filter
        x_conc = single_conc_readcount_filter(x_conc, R_filt_conc, count_cutoff)
        #calculate binding signal
        df = calc_binding_signal(
            x_conc,
            conc,
            barcodes,
            bin_key,
            weights
        )
        # merge to create full curve
        # 'outer' merge will keep all sequences even if for some
        # concentrations they have <count_cutoff
        b_curve = b_curve.merge(df, on='seq', how='outer')
    b_curve = b_curve[b_curve.isna().sum(1)<=n_missing_allowed]
    return b_curve


def get_exp_dir_dict(sample_key_file, experiment_directory):
    '''process sample_key_file to get a dictionary of experiment names
    and their paths:
    "experiment name":"path to data files"
    '''
    samples = pd.read_csv(sample_key_file)
    experiments = list(samples['experiment'].unique())
    exp_dirs = {}
    for exp in experiments:
        exp_dirs[exp] = os.path.join(experiment_directory, exp)
    return exp_dirs


def process_exp(exp_directory, experiment, sample_key, count_cutoff, excluded_barcodes, weights):
    file_list = [f for f in glob.glob(os.path.join(exp_directory, "seq_counts_complete", "*nts_only"))]
    R, T, barcode_cols = load_and_merge_data(file_list, excluded_barcodes)
    R_filtered = initial_filter(R, count_cutoff)
    read_f = read_count2fraction(R_filtered, T, barcode_cols)
    x = read_fraction2cell_count(read_f, sample_key, barcode_cols)
    b_curve = processed_data2binding_curve(x, R_filtered, weights, sample_key, experiment, count_cutoff, n_missing_allowed=2)
    return x, b_curve



def main(count_cutoff, excluded_barcodes, weights):
    experiment_directory = str(sys.argv[1])
    sample_key_file = str(sys.argv[2])

    sample_key = pd.read_csv(sample_key_file)
    exp_dirs = get_exp_dir_dict(sample_key_file, experiment_directory)
    binding_curve_folder = './binding_curves'
    if not os.path.exists(binding_curve_folder):
        os.mkdir(binding_curve_folder)
    for exp, exp_directory in exp_dirs.items():
        x, b_curve = process_exp(
            exp_directory=exp_directory,
            experiment=exp,
            sample_key=sample_key,
            count_cutoff=count_cutoff,
            excluded_barcodes=excluded_barcodes,
            weights=weights
        )
        x.to_csv('{}{}-cell_counts.csv'.format(
                os.path.join(binding_curve_folder,''),
                exp
            ), index=False)
        b_curve.to_csv('{}{}-binding_curves.csv'.format(
                os.path.join(binding_curve_folder,''),
                exp
            ), index=False)

if __name__ == "__main__":
    main(count_cutoff, excluded_barcodes, weights)
