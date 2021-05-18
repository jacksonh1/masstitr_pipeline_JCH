# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('custom_standard')
from lmfit import Model
import os
import sys



def fit_Kd(y, x):
    def basic_fit(x, init, sat, Kd):
        return init + (sat - init) * x / (x + Kd)

    gmod = Model(basic_fit)
    # print(gmod.param_names)
    # print(gmod.independent_vars)

    res = gmod.fit(list(y), x=x, init=1, sat=3000, Kd=3e-8)
    return {'Kd':res.params['Kd'].value, 'SE':res.params['Kd'].stderr}


def apply_fit(df, concentrations, cols):
    df = df.dropna()
    df_app = df[cols].apply(
        fit_Kd,
        axis=1,
        result_type='expand',
        **{'x':concentrations}
    )
    df = pd.concat([df, df_app], axis = 'columns')
    return df


def main():
    input_name = str(sys.argv[1])
    # input_name = './binding_curves/exp4-binding_curves.csv'
    sample_key_file = str(sys.argv[2])
    # sample_key_file = 'FINAL_sample_key_w_cell_counts.csv'
    output_name = os.path.join(os.path.splitext(input_name)[0], '_fit.csv')

    sample_key = pd.read_csv(sample_key_file)
    concentrations = list(sample_key['concentration (nM)'].unique())
    cols = [str(i) for i in concentrations]
    df = pd.read_csv(input_name)
    df = apply_fit(df, concentrations, cols)
    df.to_csv(output_name)

if __name__ == "__main__":
    main()

# %%
