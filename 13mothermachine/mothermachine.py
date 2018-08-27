
# coding: utf-8

# # Mother machine figure showing clpX degradation of FtsZ
#
# <a id = "toc"></a>
# ## Table of Contents
# 1. [Load data](#loaddata)
# 2. [Fluoresence ensemble plot](#ensemble)
#
# ### Load modules

# In[1]:

from __future__ import division

import numpy as np
import pandas as pd
pd.options.display.float_format = '{:,.3f}'.format
from IPython.display import display, HTML

# plotting modules and settings
import matplotlib as mpl
import matplotlib.pyplot as plt
# %matplotlib inline
import seaborn as sns
sns.set(style='ticks', color_codes=True)
sns.set_palette('deep')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Myriad Pro'

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# set axes and tick width
plt.rc('axes', linewidth=0.5)
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['xtick.minor.size'] = 2
mpl.rcParams['xtick.minor.width'] = 0.5
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['ytick.minor.size'] = 2
mpl.rcParams['ytick.minor.width'] = 0.5


# <a id = 'loaddata'></a>
# # 1. Load data
# * For each csv, each row corresponds to a single cell lineage (which undergoes growth and division over time).
# * The first two columes are the field-of-view and channel id, which uniquely identify a lineage.
# * The third and fourth columes are the absolute time index of the image and the total cellular fluorescence of the cell at that time.
# * The times are converted to be relative to the shift time using the values below.

# In[2]:

# load data frames
clpXplus_df = pd.read_csv('clpXplus.csv',sep=';')
clpXminus_df = pd.read_csv('clpXminus.csv',sep=';')

print('There are {} lineages from clpX plus data'.format(len(clpXplus_df)))
print('There are {} lineages from clpX minus data'.format(len(clpXminus_df)))

# phase picture taking interval in minutes
time_int = 2
# fluorescent picture taking interval as a function of phase images (every 30 images = 1 hr)
fl_int = 30
# absolute shift index time
clpXplus_shift_t = 63 # ~2 hours after image start
clpXminus_shift_t = 125 # ~4 hours after image start
# maximum time index
cplXplus_maxt = 500 # ~16.5 hours
cplXminus_maxt = 580 # ~19.3 hours
# normalization factor (average value from 2 hours before shift down)
clpXplus_norm = 26835
clpXminus_norm = 62009


# <a id='ensemble'></a>
# # 2. Fluorescence ensemble plot

# In[3]:

## Set up plot
fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(6, 3), squeeze=False)
ax = axes.flat[0]

# clpXplus ######################
alt_time = clpXplus_shift_t * time_int / 60 # plotting by hours realtive to shift
# arrays for holding averge values
avg_times = list(np.arange(1, cplXplus_maxt, step=fl_int))
fl_by_time = [[] for i in avg_times]

for index, row in clpXplus_df.iterrows():
    lin_times = row['lineage_times']
    # deal with formatting, removing brackets and convert to list of numbers
    lin_times = np.array([np.int(t) for t in lin_times[1:-1].split()])
    lin_fl = row['fl_totals']
    lin_fl = np.array([np.int(t) for t in lin_fl[1:-1].split()])

    # comment out to plot absolute values
    lin_fl = [x / clpXplus_norm for x in lin_fl]

    # add data to averages lists
    for i, t in enumerate(lin_times):
        idx = avg_times.index(t)
        fl_by_time[idx].append(lin_fl[i])

    # convert times and plot
    lin_times = lin_times * time_int / 60 - alt_time

    ax.plot(lin_times, lin_fl, lw=0.5, alpha=0.15, color='blue')

# average data
avg_fl = [np.mean(fl_values) for fl_values in fl_by_time]
# med_fl = [np.median(fl_values) for fl_values in fl_by_time]
n_fl = [len(fl_values) for fl_values in fl_by_time]
avg_times = np.array(avg_times) * time_int / 60 - alt_time
ax.plot(avg_times, avg_fl, color='blue', lw=1, label='clpX+ lineages average')

# clpXminus ######################
alt_time = clpXminus_shift_t * time_int / 60 # plotting by hours realtive to shift
avg_times = list(np.arange(1, cplXminus_maxt, step=fl_int))
fl_by_time = [[] for i in avg_times]

for index, row in clpXminus_df.iterrows():
    lin_times = row['lineage_times']
    # deal with formatting, removing brackets and convert to list of numbers
    lin_times = np.array([np.int(t) for t in lin_times[1:-1].split()])
    lin_fl = row['fl_totals']
    lin_fl = np.array([np.int(t) for t in lin_fl[1:-1].split()])

    # comment out to plot absolute values
    lin_fl = [x / clpXminus_norm for x in lin_fl]

    # add data to averages lists
    for i, t in enumerate(lin_times):
        idx = avg_times.index(t)
        fl_by_time[idx].append(lin_fl[i])

    # convert times and plot
    lin_times = lin_times * time_int / 60 - alt_time

    ax.plot(lin_times, lin_fl, lw=0.5, ls=':', alpha=0.15, color='green')

# average data
avg_fl = [np.mean(fl_values) for fl_values in fl_by_time]
# med_fl = [np.median(fl_values) for fl_values in fl_by_time]
n_fl = [len(fl_values) for fl_values in fl_by_time]
avg_times = np.array(avg_times) * time_int / 60 - alt_time
ax.plot(avg_times, avg_fl, color='green', lw=1, label='clpX- lineages average')

##### shared plotting and formatting
# shift down
ax.axvline(0, color='purple', ls='--', lw=0.5, label='shift-down time')

# formatting
ax.set_xlabel('time relative to shift down [hours]')
ax.set_xticks(np.arange(-2,14,2))
ax.set_xticklabels([str(x) for x in np.arange(-2,14,2)])
ax.set_xlim([-2, 10])

ax.set_ylabel('normalized fluorescence [AU]')
ax.set_yticks(np.arange(0,2.5,0.5))
ax.set_yticklabels([str(x) for x in np.arange(0,2.5,0.5)])
ax.set_ylim(bottom=0, top=1.5)

ax.legend(loc='upper right')
sns.despine()

ax.set_title('Total cellular FtsZ-mVenus fluorescence by lineage during shift down')
plt.tight_layout()
fig.savefig('./mothermachine_figure.pdf', dpi=600)
plt.show()


# In[ ]:
