# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('custom_standard')
from lmfit import Model
import os
import sys
from Bio import Seq


def trans_string(s):
    l = int(len(s)/3)*3
    return str(Seq.Seq(s[0:l]).translate())


def fit_Kd(y, x):
    def basic_fit(x, init, sat, Kd):
        return init + (sat - init) * x / (x + Kd)

    gmod = Model(basic_fit, nan_policy='omit')
    # print(gmod.param_names)
    # print(gmod.independent_vars)
    gmod.set_param_hint('Kd', value=800, min=0, max=40000)
    res = gmod.fit(y, x=x, init=1, sat=3000)
    r2 = 1 - res.redchi / np.var(y, ddof = 2)
    d = {
        'Kd':res.params['Kd'].value,
        'SE':res.params['Kd'].stderr,
        'init':res.params['init'].value,
        'sat':res.params['sat'].value,
    }
    # 'R2':r2,
    return pd.Series(d, index = ['Kd','SE','init','sat'])


def fit_long_form_data(df):
    # df = df.dropna()
    x = list(df['concentration (nM)'].astype(float))
    y = list(df['binding signal'])
    fit_params = fit_Kd(y, x)
    return fit_params


def main():
    input_name = str(sys.argv[1])
    # sample_key_file = str(sys.argv[2])
    output_name = os.path.splitext(input_name)[0] + '_fit.csv'
    # sample_key = pd.read_csv(sample_key_file)
    # concentrations = list(sample_key['concentration (nM)'].unique())
    # cols = [str(i) for i in concentrations]
    df = pd.read_csv(input_name)
    df2 = df.melt(id_vars='seq', var_name='concentration (nM)', value_name='binding signal')
    seqgrp = df2.groupby('seq')
    fits = seqgrp.apply(fit_long_form_data)
    if not fits.empty:
        fits['% fit error'] = (fits['SE']/fits['Kd'])*100
        fits = fits.reset_index()
        fin_df = pd.merge(df, fits, on='seq')
        fin_df = fin_df.sort_values(by='seq',ascending=True)
        fin_df.to_csv(output_name, index=False)
    else:
        print("binding curve dataframe was empty after trying to fit (and probably before fitting as well) most likely no sequences made it through the read/cell count filters. - check your data")


if __name__ == "__main__":
    main()

# %%
