"""
Created by Jose A. Capitan
"""

import numpy as np 
import pandas as pd
import random, math, sys, types, time
from scipy import stats, optimize
from math import erf, sqrt, fabs, pi, exp
import matplotlib.pyplot as plt
from matplotlib.mlab import find
from matplotlib import cm, rc
import statsmodels.api as sm
from glob import glob
import scipy, scipy.stats
from statsmodels.distributions.empirical_distribution import ECDF

#####    NEW STUFF FOR CLIBRARY
import ctypes
import os

permlib = ctypes.CDLL('permlib.so')

def permut_cwrap(rhovec, Smeta, Slocal, numreal):
    size=Smeta*Smeta
    sizeDoubleType = ctypes.c_double * size
    rho = sizeDoubleType()
    SlocalDoubleType = ctypes.c_double * numreal
    av_rho = SlocalDoubleType()
    
    for n in xrange(size):
        rho[n] = rhovec[n]
    
    permlib.permutations_skew(rho, ctypes.c_int(Smeta), ctypes.c_int(Slocal),
                          ctypes.c_int(numreal), av_rho)

    nprho=np.zeros(numreal)
    n=0
    for rr in av_rho:
        nprho[n]=rr
        n+=1
    return nprho

def myfind(list01):
    return np.array([n for n,a in enumerate(list01) if a], dtype=int)
    
def load_trait_matrix_ranges_stdz(fname, axis, logs=False, remove_zeros=True):
    print("Trait %d:"%(axis+1))
    # The following line should be adapted to the trait table's format; ranges of the environmental
    # factor for each species are to be placed in the last two columns
    trait = np.loadtxt(fname, delimiter=';', skiprows=1, usecols=(1,2,3,4,5,6,7,8,9,10,11,12))
    if remove_zeros:
        trait2=[]
        min2=[]
        max2=[]
        for i, row in enumerate(trait):
            if row[axis]<-9999:
                remove_sp.append(i)
            else:
                trait2.append(row[axis])
                min2.append(trait[i][10])
                max2.append(trait[i][11])
        trait=np.array(trait2)
        min=np.array(min2)
        max=np.array(max2)
    if logs:
        trait = np.log(trait)
    print("trait zeros: %d" % (sum(trait==0)))
    mu=trait.min()
    std=trait.max()-mu #Range standardization
    trait=(trait-mu)/std
    return trait,min,max

def pairwise_diffs(mat):
    df=[]
    for i,r in enumerate(mat):
        for r2 in mat:
            df.append(r-r2)
    return np.array(df)

def distance_trait(trait_diffs,nnodes):
    dist=trait_diffs
    maxdist=max(dist)
    dist/=maxdist
    maxdist=max(dist)
    distm=list2matrix(dist,nnodes)
    return(dist,distm,maxdist)

def local_traits(trait, species):
    ind = myfind(species) 
    return trait[ind]
    
def mapping(array):
    n=len(array)
    corresp=np.zeros(n)-1.
    k=0
    for i in range(n):
        if array[i]==1:
            corresp[i]=k;
            k+=1
    return(corresp)

def list2matrix(data,num_nodes):
    matrix=np.zeros((num_nodes,num_nodes))
    k=0
    for i in range(num_nodes):
        for j in range(num_nodes):
            matrix[i,j]=data[k]
            k+=1
        matrix[i,i]=0.
    return(matrix)

def submatrix(matrix,abun,corresp):
    n=len(abun)
    vec=[]
    k=0
    for i in range(n):
        if abun[i]==1:
            vec.append(corresp[i])
    vec=np.asarray(vec,dtype=int)
    k=len(vec)
    submat=np.zeros((k,k))
    for i in range(k):
        for j in range(k):
            submat[i,j]=matrix[vec[i],vec[j]]
    return(submat)

def normal_cdf(x,mu,std):
    pnorm=(1.+math.erf((x-mu)/math.sqrt(2.)/std))/2.
    return(pnorm)

def remove_negative(matrix,nnodes):
    k=0
    out=np.zeros(nnodes*(nnodes-1)/2)
    for i in range(nnodes):
        for j in range(nnodes-i-1):
            if matrix[i][j+i+1]>0.:
                out[k]=matrix[i][j+i+1]
            elif matrix[i][j+i+1]<0.:
                out[k]=matrix[j+i+1][i]
            else:
                out[k]=matrix[i][j+i+1]
            k+=1
    return(out)

def main(log_cols, num_real, axis):
    # Load trait table
    trait,min_env,max_env = load_trait_matrix_ranges_stdz('traits_plus_range.csv', axis, log_cols)
    # Load presence-absence or abundance table
    data = pd.read_csv('table_presence_absence.csv', sep='\t')
    time.sleep(1)
    # If some species have missing trait values, set their values equal to -10000 and
    # they will not be condirered in the analysis
    print("\t%s" % str(remove_sp))
    data.drop(data.index[remove_sp], inplace=True)
    # List of species
    sp_cols = list(data.index)
    # Number of species
    nspecies = len(sp_cols)
    # Output file
    if len(log_cols)>0:
        seed3 = '_log'+''.join([str(c) for c in log_cols])
    else:
        seed3 = ''
    print("\tFile Metacom_test0_ax%d%s.csv" % (axis, seed3))
    fp = open('Metacom_test0_ax%d%s.csv'%(axis, seed3),"w")
    # File headers
    ind1 = ['Dist']
    indices = [i+'.Meta,'+i+','+i+'.p' for i in ind1]
    fp.write('Sample,MetaCom.S,Local.S,'+','.join(indices)+'\n')
    # Compute cumulative abundance distribution for the species pool
    abun_meta = np.zeros(nspecies, dtype=int)
    for column in data:
        if column != 'genome':
            abun_meta = np.add(abun_meta,np.asarray(data[column]))
    for i in range(nspecies):
        if abun_meta[i] > 0:
            abun_meta[i] = 1
        else:
            abun_meta[i] = 0
    nspecies_meta = abun_meta.sum()
    # Number of species in the metacommunity
    print("\t\t%d species" % (nspecies_meta))
    corresp = mapping(abun_meta)
    # Metacommunity trait matrix
    trait_meta = local_traits(trait, abun_meta)
    # Pairwise differences
    diff_meta = pairwise_diffs(trait_meta)
    # Dissimilarity values
    dmeta_gower, matmeta_gower, mmeta_gower = distance_trait(diff_meta,nspecies_meta)
    dmeta_gower1 = remove_negative(matmeta_gower,nspecies_meta)
    # Average species pool dissimilarity
    meta_gower = dmeta_gower1.mean()
    # Loop in local communities (samples)
    for column in data:
        if column != 'genome':
            # Vector of local abundances
            abun_local = np.array(data[column], dtype=int)
            print("\t\t\t%s: %d species" % (column, abun_local.sum()))
            # p-value calculation
            matgower = submatrix(matmeta_gower, abun_local, corresp)
            dgower = remove_negative(matgower, abun_local.sum())
            # Average dissimilarity index for the sample
            gower = dgower.mean()
            dist_gower = permut_cwrap(dmeta_gower,nspecies_meta,abun_local.sum(),num_real)
            # Estimate p-value with a Gaussian cdf
            mu_gower = dist_gower.mean()
            std_gower = dist_gower.std()
            pval = normal_cdf(gower,mu_gower,std_gower)
            # Output results
            fp.write("{0},{1},{2},{3},{4},{5}\n".format(column, nspecies_meta, abun_local.sum(), meta_gower, gower, pval))
    fp.close()

# Inputs

num_real = int(sys.argv[1])   # Number of realizations of null-model distributions for p-value computation

# Optional: if some trait values are to be taken in log scale, append their indices in the last input arguments
log_cols = []
for a in sys.argv[2:]:
    log_cols.append(int(a))
print("log cols: %s" % str(log_cols))

remove_sp=[]

# Main function, loop in different traits

for axis in range(10):
    # axis -> Trait index to be analyzed (column number from 0 to n-1 in the table of traits)
    main(log_cols, num_real, axis)

#
# Usage examples from command line
#
# EXAMPLE 1: genenerate p-value distributions for all traits (0--9) using 200 null model realizations
# to calculate p-values (purely random null model)
#
# python gen-Metacom.py 200
#
# EXAMPLE 2: genenerate p-value distributions for all traits (0--9) using 100 null model realizations
# to calculate p-values (purely random null model)
#
# python gen-Metacom.py 100
#
# EXAMPLE 3: plot p-value distributions using the distributions calculated in EXAMPLE 2
#
# python boxplot png
#

