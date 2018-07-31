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

# Clibrary types
import ctypes
import os

permlib = ctypes.CDLL('permlib.so')

# Auxiliaty functions

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

def permut_cwrap_ranges(rhovec, min_env, max_env, Smeta, Slocal, env_local, numreal):
    size=Smeta*Smeta
    sizeDoubleType = ctypes.c_double * size
    rho = sizeDoubleType()
    SlocalDoubleType = ctypes.c_double * numreal
    av_rho = SlocalDoubleType()
    SmetaDoubleType = ctypes.c_double * Smeta
    min = SmetaDoubleType()
    max = SmetaDoubleType()
    
    for n in xrange(size):
        rho[n] = rhovec[n]
    
    for n in xrange(Smeta):
        min[n] = min_env[n]
        max[n] = max_env[n]

    permlib.permutations_skew_ranges(rho, min, max, ctypes.c_int(Smeta), ctypes.c_int(Slocal),
                                     ctypes.c_int(numreal), ctypes.c_double(env_local), av_rho)

    nprho=np.zeros(numreal)
    n=0
    for rr in av_rho:
        nprho[n]=rr
        n+=1
    return nprho

def permut_cwrap_weighted(rhovec, Smeta, Slocal, numreal, cumfreqvec):
    size=Smeta*Smeta
    sizeDoubleType = ctypes.c_double * size
    rho = sizeDoubleType()
    SlocalDoubleType = ctypes.c_double * numreal
    av_rho = SlocalDoubleType()
    SmetaDoubleType = ctypes.c_double * Smeta
    cumfreq = SmetaDoubleType()
    
    for n in xrange(size):
        rho[n] = rhovec[n]
    
    for n in xrange(Smeta):
        cumfreq[n] = cumfreqvec[n]
    
    permlib.permutations_skew_weighted(rho, ctypes.c_int(Smeta), ctypes.c_int(Slocal),
                                       ctypes.c_int(numreal), av_rho, cumfreq)

    nprho=np.zeros(numreal)
    n=0
    for rr in av_rho:
        nprho[n]=rr
        n+=1
    return nprho

def permut_cwrap_ranges_weighted(rhovec, min_env, max_env, Smeta, Slocal, env_local, numreal, cumfreqvec):
    size=Smeta*Smeta
    sizeDoubleType = ctypes.c_double * size
    rho = sizeDoubleType()
    SlocalDoubleType = ctypes.c_double * numreal
    av_rho = SlocalDoubleType()
    SmetaDoubleType = ctypes.c_double * Smeta
    min = SmetaDoubleType()
    max = SmetaDoubleType()
    cumfreq = SmetaDoubleType()
    
    for n in xrange(size):
        rho[n] = rhovec[n]
    
    for n in xrange(Smeta):
        min[n] = min_env[n]
        max[n] = max_env[n]
        cumfreq[n] = cumfreqvec[n]

    permlib.permutations_skew_ranges_weighted(rho, min, max, ctypes.c_int(Smeta), ctypes.c_int(Slocal),
                                              ctypes.c_int(numreal), ctypes.c_double(env_local), av_rho, cumfreq)

    nprho=np.zeros(numreal)
    n=0
    for rr in av_rho:
        nprho[n]=rr
        n+=1
    return nprho


def myfind(list01):
    return np.array([n for n,a in enumerate(list01) if a], dtype=int)
    
def load_trait_matrix_stdz(fname, axis, logs=False, remove_zeros=True):
    print("Trait %d:" % (axis+1))
    # The following line should be adapted to the trait table's format
    trait = np.loadtxt(fname, delimiter='\t', skiprows=1, usecols=(1,2,3,4,5,6,7,8,9,10))
    if remove_zeros:
        trait2=[]
        for i, row in enumerate(trait):
            if row[axis]<-9999:
                remove_sp.append(i)
            else:
                trait2.append(row[axis])
        trait=np.array(trait2)
    if logs:
        trait = np.log(trait)
    print("trait zeros: %d" % sum(trait==0))
    mu=trait.min()
    std=trait.max()-mu #Range standardization
    trait=(trait-mu)/std
    return trait

def load_trait_matrix_ranges_stdz(fname, axis, logs=False, remove_zeros=True):
    print("Trait %d:" % (axis+1))
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
    print("trait zeros: %d" % sum(trait==0))
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

def myshuffle(list):
    size = len(list)
    idx = np.random.choice(size,size,replace=False)
    shuffled_list = []
    for i in range(size):
        shuffled_list.append(list[idx[i]])
    return(shuffled_list)

def main(log_cols, num_real, axis, nrep, test, ord_flag):
    # Load trait table
    trait,min_env,max_env = load_trait_matrix_ranges_stdz('traits_plus_range.csv', axis, log_cols)
    # Load presence-absence or abundance table
    data = pd.read_csv('table_presence_absence.csv', sep='\t')
    time.sleep(1)
    # Output files
    if len(log_cols)>0:
        seed3 = '_log'+''.join([str(c) for c in log_cols])
    else:
        seed3 = ''
    if ord_flag != 1:
        ord_flag = 0
    print("File hindex_test%d_ord%d_rep%d_ax%d%s.txt" % (test, ord_flag, nrep, axis, seed3))
    fp = open('hindex_test%d_ord%d_rep%d_ax%d%s.txt'%(test, ord_flag, nrep, axis, seed3),"w")
    # File headers
    fp.write('Env.variable,Clust.Low,Clust.Med,Clust.Up,ODisp.Low,ODisp.Med,ODisp.Up\n')
    # If some species have missing trait values, set their values equal to -10000 and
    # they will not be condirered in the analysis
    print("\t%s" % str(remove_sp))
    data.drop(data.index[remove_sp], inplace=True)
    datacp = data
    # List of species
    sp_cols = list(data.index)
    # Number of species
    nspecies = len(sp_cols)
    # Sequential removal of samples in decreasing order of the environmental variable
    salinity = pd.read_csv('metadata.csv', sep='\t', usecols=[0,1], chunksize=1)
    abundances = pd.read_csv('metacomm_sum.csv', sep='\t', usecols=[1], chunksize=1)
    muestra = []
    sal = []
    for df in salinity:
        muestra.append(list(df['sample_ID']))
        sal.append(np.array(df['salinity'], dtype=float))
    nmuestras = len(muestra)
    nabun = []
    for df in abundances:
        nabun.append(np.array(df['sums'], dtype=float))
    # Initialize matrices for clustering and over-dispersion indices
    qmat = np.zeros((nrep,nmuestras))
    pmat = np.zeros((nrep,nmuestras))
    # Loop in repetitions for averages
    for rep in range(nrep):
        data = datacp
        # Random ordering of samples
        if ord_flag == 1:
            ord_muestra = myshuffle(muestra)
        else:
            ord_muestra = muestra
        print("Repetition %d" % (rep))
        # Sequential removal of samples
        for j in range(nmuestras-1):
            data = data.drop(str(ord_muestra[j][0]), axis=1)
            print("\tRemoving sample %s" % ord_muestra[j][0])
            # Compute cumulative abundance distribution for the species pool
            abun_meta = np.zeros(nspecies, dtype=int)
            for column in data:
                if column != 'genome':
                    abun_meta = np.add(abun_meta,np.asarray(data[column]))
            remove_sp1 = []
            for i in range(nspecies):
                if abun_meta[i] > 0:
                    abun_meta[i] = 1
                else:
                    abun_meta[i] = 0
                    remove_sp1.append(i)
            nspecies_meta = abun_meta.sum()
            remove_sp1 = np.array(remove_sp1)
            min_env1 = np.delete(min_env,remove_sp1,0)
            max_env1 = np.delete(max_env,remove_sp1,0)
            nabun_sum = np.zeros(nspecies_meta, dtype=float)
            i = 0
            for n in range(nspecies):
                if abun_meta[n] == 1:
                    nabun_sum[i] = nabun[n]
                    i += 1
            nabun_cum = np.zeros(nspecies_meta, dtype=float)
            nabun_cum[0] = nabun_sum[0]
            for n in range(nspecies_meta-1):
                nabun_cum[n+1] = nabun_cum[n]+nabun_sum[n+1]
            nabun_cum /= nabun_cum[nspecies_meta-1]
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
            qindex = 0
            pindex = 0
            m = 0
            for column in data:
                if column != 'genome':
                    env_local = sal[ord_muestra.index([column])][0]
                    # Vector of local abundances
                    abun_local = np.array(data[column], dtype=int)
                    # p-value calculation
                    matgower = submatrix(matmeta_gower, abun_local, corresp)
                    dgower = remove_negative(matgower, abun_local.sum())
                    # Average dissimilarity index for the sample
                    gower = dgower.mean()
                    # Randomization tests
                    if test==0: # Purely random randomization test
                        dist_gower = permut_cwrap(dmeta_gower,nspecies_meta,abun_local.sum(),num_real)
                    if test==1: # Environmental ranges randomization test
                        dist_gower = permut_cwrap_ranges(dmeta_gower,min_env1,max_env1,nspecies_meta,abun_local.sum(),env_local,num_real)
                    if test==2: # Abundace-based randomization test
                        dist_gower = permut_cwrap_weighted(dmeta_gower,nspecies_meta,abun_local.sum(),num_real,nabun_cum)
                    if test==3: # Ranges and abundances randomization test
                        dist_gower = permut_cwrap_ranges_weighted(dmeta_gower,min_env1,max_env1,nspecies_meta,abun_local.sum(),env_local,num_real,nabun_cum)
                    
                    # Estimate p-value with a Gaussian cdf
                    mu_gower = dist_gower.mean()
                    std_gower = dist_gower.std()
                    pval = normal_cdf(gower,mu_gower,std_gower)
                    
                    # Clustering and over-dispersion indices calcularion
                    if pval<0.05:
                        qindex += 1
                    if pval>0.95:
                        pindex += 1
                    m += 1
            print("\th-index: %lf" % float(float(qindex)/float(m)))
            qmat[rep][j] = float(qindex)/float(m)
            pmat[rep][j] = float(pindex)/float(m)

    # Output results
    sal_tot = 0.
    for j in range(nmuestras):
        sal_tot += sal[j][0]
    for j in range(nmuestras-1):
        qvec = qmat[:,j]
        pvec = pmat[:,j]
        # Clustering (p < 0.05) and overdispersion indices (p > 0.95) are stored
        fp.write("{0},{1},{2},{3},{4},{5},{6}\n".format(sal_tot/float(nmuestras-j+1), np.percentile(qvec,5), np.percentile(qvec,50), np.percentile(qvec,95), np.percentile(pvec,5), np.percentile(pvec,50), np.percentile(pvec,95)))
        fp.flush()
        sal_tot -= sal[j][0]

    fp.close()

# Inputs

num_real = int(sys.argv[1])   # Number of realizations of null-model distributions for p-value computation
axis = int(sys.argv[2])       # Trait index to be analyzed (column number from 0 to n-1 in the table of traits)
nrep = int(sys.argv[3])       # Number of repetitions to average h index
test = int(sys.argv[4])       # Randomization test indicatior (0,1,2, or 3)
ord_flag = int(sys.argv[5])   # Random ordering (1) or decreasing salinity order (other value)

# Optional: if some trait values are to be taken in log scale, append their indices in the last input arguments
log_cols = []
for a in sys.argv[6:]:
    log_cols.append(int(a))
print("log cols: %s" % str(log_cols))

remove_sp = []

# Main function
main(log_cols, num_real, axis, nrep, test, ord_flag)

#
# Usage examples from command line
#
# EXAMPLE 1: genenerate h index for trait 3 (0--9) using 100 realizations for null model 2 (weighted by abundances)
# with decreasing order of the environmental variable (ord_flag = 0). h index is averaged using 50 repetitions
#
# python hindex-repetitions.py 100 3 50 2 0
#
# EXAMPLE 2: genenerate h index for trait 7 (0--9) using 200 realizations for null model 1 (weighted by ranges of the
# environmental gradient) with decreasing order of the environmental variable (ord_flag = 0). h index is averaged
# using 50 repetitions
#
# python hindex-repetitions.py 200 7 50 1 0
#
# EXAMPLE 3: genenerate h index for trait 1 (0--9) using 100 realizations for null model 3 (weighted by anbundances
# and ranges of the environmental gradient) using random orderings of the environmental variable (ord_flag = 1).
# h index is averaged over 100 repetitions
#
# python hindex-repetitions.py 100 1 100 3 1
#


