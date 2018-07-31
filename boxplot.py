"""
Created by Jose A. Capitan
"""

import numpy as np
import random, math, sys, types
from scipy import stats, optimize
from math import erf, sqrt, fabs, pi, exp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.mlab import find
from matplotlib import cm, rc, font_manager
import statsmodels.api as sm
import os
from glob import glob
import scipy, scipy.stats
from scipy.stats.stats import pearsonr
from matplotlib.patches import Polygon
import matplotlib.ticker as ticker

## PLOT DEFAULTS
################

rc('xtick', labelsize=40, color='k')
rc('ytick', labelsize=40, color='k')
rc('text', fontsize=50)
rc('axes', labelsize=40)
rc('axes', edgecolor='k')
rc('axes', titlesize=50)
rc('legend', fontsize='x-large')
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{sfmath}')

def iterable(obj):
    return hasattr(obj, '__iter__')

def load_data(f='Metacom_log23.Woody.csv'):
    return pd.read_csv(f,sep=',')

def get_group_id(gdata, i):
    n=1
    for name, group in gdata:
        if (n==i):
            return name, group
        n+=1

def get_pval(group, index, metaindex):
    x=np.array(group[index])
    y=group[metaindex].mean()
    S=group['MetaCom.S'].mean()
    return x, y, S

def gen_box_colors(size, cmap=cm.Paired):
    return cmap(np.linspace(0, 1, size))

def write_box_colors(path, names, carray):
    fp=open(path, "w")
    cols = np.array(carray*256, dtype=int)
    for n, (name, col) in enumerate(zip(names, cols)):
        fp.write('{0} {4} ({1}, {2}, {3})\n'.format(n, col[0], col[1], col[2], name))
    fp.close()

def write_box_colors2(path, carray):
    fp=open(path, "w")
    cols = np.array(carray*256, dtype=int)
    for n, col in enumerate(cols):
        #@map color 0 to (255, 255, 255), "white"
        fp.write('@map color {0} to ({1},{2},{3}), "{4}"\n'.format(20+n, col[0], col[1], col[2], n))
    fp.close()

def boxplot_by_trait(filename, filext):

    names=['Genome Size','$\%$ GC','$\%$ Coding base','$\%$ CDS','$\%$ RNA','rRNA count','$\%$ Transporter','$\%$ Signal Peptide','$\%$ Transmembrane','Gene Count']
    
    data=[]
    n=1
    for trait in range(10):
        csvdata=load_data('%s_ax%d.csv'%(filename,trait))
        data.append(csvdata['Dist.p'])
        n+=1
    
    fig, ax1 = plt.subplots(figsize=(20,15))
    plt.subplots_adjust(left=0.10, right=0.985, top=0.96, bottom=0.35)

    bp = plt.boxplot(data, notch=0, sym='k.', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black', linewidth=2)
    plt.setp(bp['whiskers'], color='black', linestyle='-', linewidth=2)
    plt.setp(bp['caps'], color='black', linestyle='-', linewidth=2)
    ax1.set_axisbelow(True)
    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    ax1.yaxis.set_major_locator(yloc)
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    xtickNames = plt.setp(ax1, xticklabels=names)
    plt.setp(xtickNames, rotation=90)
    ax1.tick_params('both', length=10)
    ax1.set_ylabel('p-value distribution')
    boxColors = ['royalblue','white','lightgray']

    spanColors = ['grey', 'white']
    numBoxes = n-1
    boxColors2 = gen_box_colors(numBoxes)
    x=0.5
    flag=1
    alpha=0.3
    plt.axhspan(0.95, 0.999999, facecolor='tomato', alpha=alpha, linewidth=0)
    plt.axhspan(0, 0.05, facecolor='steelblue', alpha=alpha, linewidth=0)
    for i in range(numBoxes):
        plt.axvspan(x, x+1, facecolor=spanColors[flag], alpha=alpha, linewidth=0)
        x+=1
        flag=1-flag
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = zip(boxX,boxY)
        boxPolygon = Polygon(boxCoords, facecolor=boxColors2[i])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k', linewidth=2)

    plt.ylim(0,1)

    plt.savefig('./boxplot_bytrait.%s'%(filext))

ext = sys.argv[1]

boxplot_by_trait('Metacom_test0', ext)

# EXAMPLE 3: plot p-value distributions using the distributions calculated in EXAMPLE 2
#
# python boxplot png
#