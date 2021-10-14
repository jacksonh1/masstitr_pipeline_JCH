import sys
import os
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use("custom_standard")


# %%


def filter_mT_table(df, kd_up_lim, SE_upper_lim, kd_low_lim=0, drop_dup=True):
    """
    filters existing masstitr table
    filters out seqs containing * or X characters
    removes seqs not matching xxxPxExxx motif

    Parameters
    ----------
    SE_upper_lim
        Will filter data for SE values with SE <= `SE_upper_lim`
    kd_up_lim
        Will remove data outside range: Kd < `kd_up_lim`
    kd_low_lim
        Will remove data outside range: `kd_low_lim` < Kd
    drop_dup : :obj:`bool`, optional
        Whether or not to drop duplicate sequences

    Returns
    ------
    filtered dataframe
    """
    # drop duplicates first just in case only a single duplicate makes the cutoff but the other(s) do not
    # for example, if LNLPEESDW had 2 entries with diff Kd values with one of them > kd_up_lim and the other < kd_up_lim
    if drop_dup:
        df = df.drop_duplicates(subset="AA_seq", keep=False)
    df = df[~df.AA_seq.str.contains("[\*X]")]
    df = df[df.AA_seq.str.contains(r"...P.E...")]
    df = df[df["SE"] <= SE_upper_lim]
    # df = df[df["R2"] >= R2_lower_lim]
    df = df[(df["Kd"] < kd_up_lim) & (df["Kd"] > kd_low_lim)]
    return df




def correlation(df1, df2, plots=False):
    """
    correlation plot and calculation between masstitr experiments
    correltaion between columns labelled `Kd`.

    Parameters
    ----------
    df1
        first dataframe
    df2
        second dataframe

    Returns
    -------
    p_r
        pearson correlation R value
    s_r
        spearman correlation R value
    n_points
        number of points in correlation
    """
    df = pd.merge(df1, df2, on='AA_seq', suffixes=('_1','_2'))
    df = df.reset_index(drop=True)
    if plots:
        fig, ax = plt.subplots(figsize=[4, 4])
        df.plot.scatter(x='Kd_1', y='Kd_2', loglog=True, ax=ax)
        fig, ax = plt.subplots(figsize=[4, 4])
        df.plot.scatter(x='Kd_1', y='Kd_2', ax=ax)
    p_r = df.corr("pearson").loc['Kd_1','Kd_2']
    s_r = df.corr("spearman").loc['Kd_1','Kd_2']
    n_points = len(df)
    return p_r, s_r, n_points


def correlation_screen(df1, df2, p_fit_err_cutoffs, kd_cutoffs):
    p_mat = pd.DataFrame(columns=p_fit_err_cutoffs, index=kd_cutoffs)
    s_mat = pd.DataFrame(columns=p_fit_err_cutoffs, index=kd_cutoffs)
    n_mat = pd.DataFrame(columns=p_fit_err_cutoffs, index=kd_cutoffs)
    for p_fit_err in p_fit_err_cutoffs:
        for Kd in kd_cutoffs:
            df1_f = filter_mT_table(
                df1,
                kd_up_lim=Kd,
                SE_upper_lim=Kd,
                # R2_lower_lim=0,
                drop_dup=True
            )
            df1_f = df1_f[df1_f['% fit error']<=p_fit_err]
            df2_f = filter_mT_table(
                df2,
                kd_up_lim=Kd,
                SE_upper_lim=Kd,
                # R2_lower_lim=0,
                drop_dup=True
            )
            df2_f = df2_f[df2_f['% fit error']<=p_fit_err]
            p, s, n = correlation(df1_f, df2_f, plots=False)
            p_mat.loc[Kd, p_fit_err] = p
            s_mat.loc[Kd, p_fit_err] = s
            n_mat.loc[Kd, p_fit_err] = n
    return p_mat, s_mat, n_mat




# %%
# ==============================================================================
# // TITLE
# ==============================================================================

def driver(directory):
    exp1_file = os.path.join(directory,"exp1-binding_curves_fit.csv")
    exp2_file = os.path.join(directory,"exp2-binding_curves_fit.csv")
    df1 = pd.read_csv(exp1_file)
    df1 = df1.drop(['3000', '1000', '300', '100', '30', '10', '3', '1'], axis=1)
    df2 = pd.read_csv(exp2_file)
    df2 = df2.drop(['3000', '1000', '300', '100', '30', '10', '3', '1'], axis=1)


    p_fit_err_cutoffs = np.linspace(15,50,10)
    p_fit_err_cutoffs = [round(i) for i in p_fit_err_cutoffs]

    kd_cutoff = 3000
    kd_cutoffs = np.linspace(300,3000,10)
    kd_cutoffs = [round(i) for i in kd_cutoffs]

    p_mat, s_mat, n_mat = correlation_screen(df1, df2, p_fit_err_cutoffs, kd_cutoffs)

    fig, ax = plt.subplots(figsize=[9,8])
    sns.heatmap(p_mat.astype(float),vmin=0, vmax=1, annot=True, fmt=".2g", linewidth=0.5, ax=ax, annot_kws={"size": 10})
    ax.tick_params(axis=u'both', which=u'both',length=0,labelsize=16)
    ax.set_xlabel('percent error cutoff')
    ax.set_ylabel('Kd cutoff')
    ax.set_title('pearson R')
    plt.tight_layout()
    plt.savefig(os.path.join(directory,'pearson_heatmap.jpg'),format='jpg',dpi=300)

    fig, ax = plt.subplots(figsize=[9,8])
    sns.heatmap(s_mat.astype(float),vmin=0, vmax=1, annot=True, fmt=".2g", linewidth=0.5, ax=ax, annot_kws={"size": 10})
    ax.tick_params(axis=u'both', which=u'both',length=0,labelsize=16)
    ax.set_xlabel('percent error cutoff')
    ax.set_ylabel('Kd cutoff')
    ax.set_title('spearman R')
    plt.tight_layout()
    plt.savefig(os.path.join(directory,'spearman_heatmap.jpg'),format='jpg',dpi=300)

    fig, ax = plt.subplots(figsize=[9,8])
    sns.heatmap(n_mat.astype(float), annot=True, fmt=".4g", linewidth=0.5, ax=ax, annot_kws={"size": 10})
    ax.tick_params(axis=u'both', which=u'both', length=0, labelsize=16)
    ax.set_xlabel('percent error cutoff')
    ax.set_ylabel('Kd cutoff')
    ax.set_title('n points in correlation')
    plt.tight_layout()
    plt.savefig(os.path.join(directory,'n_points_heatmap.jpg'),format='jpg',dpi=300)

# %%

# %%

def main():
    directory = str(sys.argv[1])
    driver(directory)


if __name__ == "__main__":
    main()
