import re
from glob import glob
import pandas as pd
import numpy as np
from matplotlib import pyplot as pl
import seaborn as sns
from matplotlib import rcParams

# dictionary to init dataframe that stores all the information
df = {}
df['toolkit'] = []
df['n cells'] = []
df['step'] = []
df['total memory (GB)'] = []
df['memory_change'] = []
df['CPU time (min)'] = []

# match html by R output
for filename in glob('*cellranger_R*.html'):
    n_cells = int(re.findall('[0-9.]+', filename)[2])
    textfile = open(filename, 'r')
    filetext = textfile.read()
    textfile.close()
    # memory
    matches = re.findall('[0-9.]+ .B<', filetext)
    mem_tot = []
    mem_change = [0]
    for im, m in enumerate(matches):
        num = float(re.findall('[0-9.]+', m)[0])
        if 'MB' in m: num /= 1000
        if im == 0: mem_tot += [num]
        else:
            if im % 2 == 1:
                mem_tot += [num]
            else: mem_change += [num]
    df['total memory (GB)'] += mem_tot
    df['memory_change'] += mem_change
    # cpu time
    matches = re.findall('elapsed:.+ [0-9.]+', filetext)
    df['CPU time (min)'] += [0]
    df['CPU time (min)'] += [float(re.findall('[0-9.]+', m)[0])/60 for m in matches]
    # update step and all other fields
    df['memory_change'] += [np.nan, np.nan]
    df['total memory (GB)'] += [np.nan, np.nan]
    df['CPU time (min)'] += [np.nan, np.nan]
    # type
    df['toolkit'] += ['Cell Ranger (R)' for i in range(7)]
    # general
    df['step'] += ['init', 'load', 'Preprocessing', 'PCA', 'tSNE', 'diffmap', 'DPT']
    df['n cells'] += [n_cells for i in range(7)]


def extract_minutes(string):
    match_min = re.findall('[0-9.]+ *min', string)
    match_s = re.findall('[0-9.]+ *s', string)
    match_ms = re.findall('[0-9.]+ *ms', string)
    min = float(match_min[0].replace('min', '').strip()) if match_min else 0
    min += float(match_s[0].replace('s', '').strip())/60 if match_s else 0
    min += float(match_ms[0].replace('ms', '').strip())/60/1000 if match_ms else 0
    return min

def extract_memory(string):
    l = string.split(' GB, difference ')
    return float(l[0]), float(l[1].replace(' GB', ''))

# match html by Scanpy output
for filename in glob('*cellranger_Py*.html'):
    n_cells = int(re.findall('[0-9.]+', filename)[2])
    textfile = open(filename, 'r')
    filetext = textfile.read()
    textfile.close()
    # cpu time
    matches = re.findall('Wall time:.+', filetext)
    df['CPU time (min)'] += [0]
    df['CPU time (min)'] += [extract_minutes(m) for m in matches]
    # memory
    matches = re.findall('[0-9.]+ GB, difference [+\-0-9.]+ GB', filetext)
    df['total memory (GB)'] += [extract_memory(m)[0] for m in matches]
    df['memory_change'] += [extract_memory(m)[1] for m in matches]
    # type
    df['toolkit'] += ['Scanpy (Py)' for i in range(7)]
    # general
    df['step'] += ['init', 'load', 'Preprocessing', 'PCA', 'tSNE', 'diffmap', 'DPT']
    df['n cells'] += [n_cells for i in range(7)]


df = pd.DataFrame(df)
df_single = df.loc[df['toolkit'] != 'Cell Ranger (R)']
df_single = df_single.loc[df['step'] != 'init']
df_single = df_single.loc[df['step'] != 'load']
df_single = df_single.loc[df['step'] != 'Preprocessing']
df_single = df_single.loc[df['step'] != 'PCA']
df_single = df_single.loc[df['step'] != 'tSNE']

# remove uninteresting steps
df = df.loc[df['step'] != 'init']
df = df.loc[df['step'] != 'load']
df = df.loc[df['step'] != 'DPT']
df = df.loc[df['step'] != 'diffmap']

g = sns.FacetGrid(df, col='step', hue='toolkit', sharey=False, legend_out=True)
g = g.map(pl.scatter, 'n cells', 'total memory (GB)')
pl.subplots_adjust(top=0.82, right=0.82)
pl.legend(bbox_to_anchor=(1.04, 0.5), loc=2, borderaxespad=0.)
g.fig.suptitle('Process memory after step')
pl.savefig('figs/memory.png', dpi=400)

# g = sns.FacetGrid(df, col='step', hue='toolkit', sharey=False)
# g = g.map(pl.scatter, 'n cells', 'memory_change')
# pl.subplots_adjust(top=0.82, right=0.82)
# pl.legend(bbox_to_anchor=(1.04, 0.5), loc=2, borderaxespad=0.)
# g.fig.suptitle('changed memory during step (GB)')

g = sns.FacetGrid(df, col='step', hue='toolkit', sharey=False)
g = g.map(pl.scatter, 'n cells', 'CPU time (min)')
pl.subplots_adjust(top=0.82, right=0.82)
pl.legend(bbox_to_anchor=(1.04, 0.5), loc=2, borderaxespad=0.)
g.fig.suptitle('CPU time of step')
pl.savefig('figs/cpu_time.png', dpi=400)

# compute Speedup
df1 = pd.DataFrame()
df1['step'] = df['step'][df['toolkit'] == 'Scanpy (Py)'].values
df1['n cells'] = df['n cells'][df['toolkit'] == 'Scanpy (Py)'].values
df1['Speedup'] = df['CPU time (min)'][df['toolkit'] == 'Cell Ranger (R)'].values / df['CPU time (min)'][df['toolkit'] == 'Scanpy (Py)'].values
df1['memory ratio'] = 1 / df['total memory (GB)'][df['toolkit'] == 'Cell Ranger (R)'].values * df['total memory (GB)'][df['toolkit'] == 'Scanpy (Py)'].values

g = sns.FacetGrid(df1, col='step', sharey=False)
g = g.map(pl.scatter, 'n cells', 'Speedup', color='grey')
pl.subplots_adjust(top=0.82, right=0.82)
g.fig.suptitle('Speedup Scanpy vs. Cell Ranger (Zheng el al., 2017)')
pl.savefig('figs/speedup.png', dpi=400)
pl.savefig('figs/speedup.pdf')

g = sns.FacetGrid(df1, col='step', sharey=False)
g = g.map(pl.scatter, 'n cells', 'memory ratio', color='grey')
pl.subplots_adjust(top=0.82, right=0.82)
g.fig.suptitle('Memory ratio Scanpy vs. Cell Ranger (Zheng el al., 2017)')
pl.savefig('figs/memory_ratio.png', dpi=400)

# scaling DPT and diffmap
g = sns.FacetGrid(df_single, col='step', sharey=False)
g = g.map(pl.scatter, 'n cells', 'CPU time (min)')
pl.subplots_adjust(top=0.82, right=0.82)
g.fig.suptitle('CPU time Diffmap and DPT')
pl.savefig('figs/cpu_time_dpt.png', dpi=400)
