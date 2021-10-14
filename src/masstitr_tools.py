import pandas as pd
import os


def get_exp_dir_dict(sample_key_file, experiment_directory):
    '''process sample_key_file to get a dictionary of experiment names
    and the paths to their associated and sorted barcode files:
    "experiment name":"path to data files"
    '''
    samples = pd.read_csv(sample_key_file)
    experiments = list(samples['experiment'].unique())
    exp_dirs = {}
    for exp in experiments:
        exp_dirs[exp] = os.path.join(experiment_directory, exp)
    return exp_dirs


# %%
# ==============================================================
# // import data 
# ==============================================================

# TODO: add documentation
def df_import_1(file):
    '''
    imports barcode_x_f_nts_only file.
    returns a dataframe of the nt sequences (`seq`) and their counts (column named after file name `barcode_x`)
    '''
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
    DataFrame will be the nt sequences (`seq`) and their counts for each barcode file in `file_list`
    each barcode file in `file_list` should correspond to 1 column in the Dataframe. columns are sorted
    in numerical order (ex: `seq` | `barcode_1` | `barcode_2` | ...)
    """
    # R will be the read counts for each sequence in each gate
    R = pd.DataFrame(columns=['seq'])
    for f in file_list:
        if os.path.splitext(os.path.basename(f))[0].replace('_f_nts_only', '') in excluded_barcodes:
            continue
        print("processing file: {}".format(f))
        df1 = df_import_1(f)
        R = pd.merge(R, df1, on='seq', how='outer')
    R = R.fillna(0)
    # reorder columns
    barcode_cols = [col for col in R.columns]
    barcode_cols.sort()
    barcode_cols.remove('seq')
    R = R[['seq'] + barcode_cols]
    # `T` will be the total read counts in each barcode_x sample (i.e. total reads in each gate)
    T = R[barcode_cols].sum()
    return R, T, barcode_cols


# %%
# ==============================================================
# // read count normalization 
# ==============================================================


# normalization --------------------------------------------------
def read_count2fraction(R, T, barcode_cols):
    '''
    normalize by dividing read counts for each sequence by total number of
    reads in the sample.
    '''
    # read fraction: R/T
    read_f = R.copy() # probably unnecessary but let's keep the original dataframe anyway
    for col in barcode_cols:
        read_f[col] = R[col]/T[col]
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

# %%
# ==============================================================
# // read count filtering and read clustering
# ==============================================================


def trans_string(s):
    l = int(len(s)/3)*3
    return str(Seq.Seq(s[0:l]).translate())


def cluster_AA_seqs(R, barcode_cols):
    '''
    translates the nt sequences in `R` and then sums the read counts in each gate (barcode column)
    for each amino acid sequence. S
    '''
    clustered_R=R.groupby('AA_seq')[barcode_cols].sum()
    clustered_R = clustered_R.reset_index()
    clustered_R = clustered_R.rename(columns={'AA_seq':'seq'})
    clustered_R = clustered_R[['seq'] + barcode_cols]
    clustered_R['AA_seq'] = clustered_R['seq']
    return clustered_R


def initial_filter(R, count_cutoff):
    '''
    NOT USED IN FINAL PIPELINE
    filter out any sequences where all of the bins are below cutoff threshold
    i.e. keep any sequences where at least one of the barcodes has >100 reads
    only use this if `filter_method='reads'`
    '''
    return R[(R.drop("seq", axis=1) > count_cutoff).any(1)]


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


def single_conc_cell_count_filter(x_df, count_cutoff):
    """filter out sequences that don't have at least `count_cutoff` cells in one of the barcodes of `x_df`

    Parameters
    ----------
    x_df : pandas DataFrame
        DataFrame of the normalized/processed data (x) corresponding to
        the gates in a single concentration point. `x_df` should also
        have a column of the sequences labelled "seq".
        This is usually a subset of the full x DataFrame
        example: `x[['seqs', 'barcode_168', 'barcode_169', ... ]]
    count_cutoff : int
        The cell count cutoff to filter the sequences by

    Returns
    -------
    pandas DataFrame
        dataframe containing only the sequences passing the filter
    """
    return x_df[(x_df.drop("seq", axis=1) >= count_cutoff).any(1)]



# %%
# ==============================================================
# // group barcodes into gates for each concentration
# // Calculate binding signal for each concentration
# ==============================================================


# barcode organization --------------------------------------------- 
# using the sample key to assign the correct bins, exp, and concentrations to the correct barcodes
# it all gets used when calculating the binding signal
def conc2barcodes(sample_key, exp, conc):
    """
    uses the `sample_key` to get the barcodes (returned list of barcodes) corresponding to a given `exp` and `conc`
    """
    return list(sample_key[(sample_key['experiment']==exp) & (sample_key['concentration (nM)']==conc)]['barcode'])


def barcodes2bin(sample_key, barcodes):
    """
        gets the gate number corresponding to each barcode_x (found in the sample key) and generates a dictionary (`bin_key`) that is used to assign the read counts to the correct gates when calculating binding signal
        """
    bin_key = {}
    for barcode in barcodes:
        b = sample_key[sample_key['barcode']==barcode]['bin'].values[0]
        bin_key[b] = barcode
    return bin_key
