#!/usr/bin/env python3
# This script is included in documentation. Adapt line numbers if touched.

import glob
import pandas, seaborn

results = glob.glob("result_mpi_*.csv")

df_list = []
for result in results:
    df_list.append(pandas.read_csv(result))

df = pandas.concat(df_list)
seaborn.relplot(data=df, kind="line", x="t/ms", y="U/mV",hue="Cell",ci=None).savefig('mpi_result.svg')