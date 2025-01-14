#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import tabulate
from matplotlib.ticker import LogLocator
from matplotlib.colors import to_rgba
import re

# Limit matplot lib print precision
pd.set_option("display.precision", 2)
pd.set_option("display.float_format", "{:0.2f}".format)

Q = 10000

df = pd.read_json("stats.json")
df["seq_sz"] = df["seq_size_MB"]
df["sketch_sz"] = df["sketch_size_MB"]
df["index_sz"] = df["index_size_MB"]
df["total_sz"] = df["total_size_MB"]
df["mism/q"] = (df["query_mismatches"] // df["queries"]).fillna(0).astype("Int64")
df["remap"] = df["sketch_remap"].fillna(-1).astype(int)
df["k"] = df["sketch_k"].fillna(-1).astype(int)

# FIX k (TODO drop)
df.loc[df.k == 32, "k"] = 4
df.loc[df.k == 64, "k"] = 8
df.loc[df.k == 128, "k"] = 15
df.loc[df.k == 256, "k"] = 28

# Remap=None becomes Remap=0 and Remap=1
# (but not for SIndex)
dfr = df[(df.remap == -1) & (df["index"] != "SIndex")]
dfr["remap"] = 0
df = pd.concat([df, dfr])
df.loc[df.remap == -1, "remap"] = 1

# Squash two SDSL FM indices together
df.loc[df["index"] == "sdsl_lite_fm::FmIndexByte32Ptr", "index"] = "SDSL FM"
df.loc[df["index"] == "sdsl_lite_fm::FmIndexInt32Ptr", "index"] = "SDSL FM"

# Drop SDSL without remap
df = df[~((df["index"] == "SDSL FM") & (df["remap"] == 0))]


df["l"] = df["sketch_l"].fillna(-1).astype(int)
df["store_ms"] = df["index_store_ms_seq"].fillna(-1).astype(int)
df["sa_sampl"] = df["index_sa_sampling"].fillna(-1).astype(int)
df["int_width"] = df["index_width"].fillna(-1).astype(int)

df["sketch%"] = (df["t_query_sketch"] / df["query_time"] * 100).astype(int)
df["search%"] = (df["t_query_search"] / df["query_time"] * 100).astype(int)
df["check%"] = (df["t_query_check"] / df["query_time"] * 100).astype(int)
df["invert%"] = (df["t_query_invert_pos"] / df["query_time"] * 100).astype(int)
df["us/q"] = df["query_time"] / Q * 10**6
df["search/q"] = df["t_query_search"] / Q * 10**6
df["Total"] = df["Build"] + df["Sketch"]
# print(df.columns)

# fmt:off
cols = [
    # "sketch_params", "index_params",
    'index',
    'k', 'l', 'remap',
    'store_ms', 'sa_sampl', 'int_width',
    'sketch_sz', 'index_sz', 'total_sz',
    #
    'Sketch', 'Build',
    #
    'us/q', 'sketch%', 'search%', 'check%', 'invert%',
    'mism/q',
]

def kl(row):
    if row.k==-1:
        return 'Plain text index'
    else:
        return f"U-index (k,l) = ({row.k}, {row.l})"
df['kl'] = df.apply(kl, axis=1)

def params(row):
    idx = row['index']
    if row['index'] == 'SIndex':
        return idx
    r = '' if row['remap'] == 1 else ' -H'
    ms = '' if row['store_ms'] == 1 else ' -S'
    if row['store_ms'] == -1:
        return f"{idx}{r}"
    else:
        return f"{idx}{r}{ms}"
df['params'] = df.apply(params, axis=1)

order = [
    'libsais',
    'libsais -H',
    'libsais -S',
    'libsais -H -S',
    'SDSL FM',
    'AWRY',
    'AWRY -H',
    'SIndex',
]
df['order'] = df['params'].apply(lambda x: order.index(x))

df = df.groupby(['order', 'k']).first().reset_index()
df.sort_values(by=['order', 'k'], inplace=True)

dfc = df[cols]

# print(dfc)
# print(tabulate.tabulate(dfc, headers=dfc.columns, tablefmt="orgtbl", floatfmt=".1f"))
# print('Seq size: ', df['seq_sz'].unique())
# print('Matches', df['query_matches'].unique())
# print('#minimizers', df[['sketch_k', 'sketch_l', 'num_minimizers']].drop_duplicates().dropna())

df['test'] = df['total_sz'] / 2

# 3x3 grid of figures
fig, axs = plt.subplots(3, 1, figsize=(13, 11))
g1 = sns.barplot(data=df, hue='kl', y='total_sz', x='params', ax=axs[0], legend=False)
gx = sns.barplot(data=df, hue='kl', y='sketch_sz', x='params', ax=axs[0], alpha=0.3, color='black', legend=False)
g2 = sns.barplot(data=df, hue='kl', y='Total', x='params', ax=axs[1], legend=False)
gy = sns.barplot(data=df, hue='kl', y='Sketch', x='params', ax=axs[1], alpha=0.3, color='black', legend=False)
g3 = sns.barplot(data=df, hue='kl', y='us/q', x='params', ax=axs[2])
gz = sns.barplot(data=df, hue='kl', y='search/q', x='params', ax=axs[2], alpha=0.3, color='black', legend=False)
# g4 = sns.barplot(data=df, hue='kl', y='int_width', x='params', ax=axs[3])

legend = g3.legend(loc='lower center', bbox_to_anchor=(.5,-0.30), ncol=5, title='', frameon=False)

# Add one black box to the legend
g1.legend([plt.Rectangle((0,0),1,1,fc="black", alpha=0.5, edgecolor = 'none')], ['Size of minimizer positions and remap'], loc='upper right')
g2.legend([plt.Rectangle((0,0),1,1,fc="black", alpha=0.5, edgecolor = 'none')], ['Time sketching the input'], loc='upper right')
g3.legend([plt.Rectangle((0,0),1,1,fc="black", alpha=0.5, edgecolor = 'none')], ['Time spent in inner LOCATE'])
fig.add_artist(legend)

gx.set_ylim(2**1, 2**11)
gy.set_ylim(2**-2, 2**7)
gz.set_ylim(2**0, 2**16)
g1.set_yscale("log", base=2)
g2.set_yscale("log", base=2)
g3.set_yscale("log", base=2)

# Minor ticks at every power of 2
g1.yaxis.set_minor_locator(LogLocator(base=2.0, subs='all',numticks=20))
g2.yaxis.set_minor_locator(LogLocator(base=2.0, subs='all',numticks=20))
g3.yaxis.set_minor_locator(LogLocator(base=2.0, subs='all',numticks=20))

# horizontal grid on
g1.grid(axis='y', linewidth=0.5)
g2.grid(axis='y', linewidth=0.5)
g3.grid(axis='y', linewidth=0.5)
# g4.grid(axis='y', linewidth=0.5)
g1.grid(axis='y', linewidth=0.2, which='minor')
g2.grid(axis='y', linewidth=0.2, which='minor')
g3.grid(axis='y', linewidth=0.2, which='minor')


g1.set_xlabel('')
g2.set_xlabel('')
g3.set_xlabel('')
# g3.set_xlabel('Index.   -H: skip remap, -S: implicit minimizer seq')
g1.set_ylabel('Size (MB)')
g2.set_ylabel('Build (s)')
g3.set_ylabel('Query (us)')


fig.savefig(f"plot.svg", bbox_inches="tight")
# fig.savefig(f"plot.png", bbox_inches="tight", dpi=400)
# plt.show()
