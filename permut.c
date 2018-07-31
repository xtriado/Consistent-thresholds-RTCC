#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Commands to compile the library permlib.so
// gcc -shared -o permlib.so -fPIC permut.c -lc -lm

double get_random() { return ((double)rand() / (double)RAND_MAX); }

void sample_norep(int size, int maxi, int *vector);
int lookfor(int *vector, int value, int size);
void sample_norep_weighted(int size, int maxi, int *vector, double *cumfreq);
int sample_weighted(double unif, int size, double *cumfreq);
void permutations(double *rhovec, int Smeta, int Slocal, int num_real, double *av_rho);
void permutations_skew(double *rhovec, int Smeta, int Slocal, int num_real, double *av_rho);
void permutations_skew_weighted(double *rhovec, int Smeta, int Slocal, int num_real, double *av_rho,
                                double *cumfreq);
void permutations_weighted(double *trait, int Smeta, int Slocal, int num_real, double *av_trait, double *cumfreq);
void permutations_skew_ranges(double *rhovec, double *min_env, double *max_env, int Smeta, int Slocal, double env_local, int num_real, double *av_rho);
void permutations_skew_ranges_weighted(double *rhovec, double *min_env, double *max_env, int Smeta, int Slocal, double env_local, int num_real, double *av_rho, double *cumfreq);
int sample_norep_ranges(double *min_env, double *max_env, int size, int maxi, double env_local, int *vector);
int sample_norep_ranges_weighted(double *min_env, double *max_env, int size, int maxi, double env_local, int *vector, double *cumfreq);

void permutations_weighted(double *trait, int Smeta, int Slocal, int num_real, double *av_trait, double *cumfreq)
{
    int i, j, *local_com_indices = malloc(Slocal*sizeof(int));
    double rval;
    time_t tp;
    unsigned long cputime;

    time(&tp);
    cputime = (unsigned long)tp;

    srand(cputime);
    
    for (i = 0; i < num_real; i++) {
        
        sample_norep_weighted(Slocal, Smeta, local_com_indices, cumfreq);
        
        av_trait[i] = 0.;
        
        for (j = 0; j < Slocal; j++) {
            
            rval = trait[local_com_indices[j]];
            av_trait[i] += rval;
        }
        
        av_trait[i] /= (double)(Slocal);
    }
    
    free(local_com_indices);
    
    return;
}

void permutations(double *rhovec, int Smeta, int Slocal, int num_real, double *av_rho)
{
    int i, j, k, *local_com_indices = malloc(Slocal*sizeof(int));
    // rho must be the full symmetric metacommunity matrix, expressed as a vector with
    // diagonal elements equal to one
    double **rho;
    time_t tp;
    unsigned long cputime;
    
    time(&tp);
    cputime = (unsigned long)tp;
    
    srand(cputime);
    
    k = 0;
    rho = (double **) malloc(Smeta*sizeof(double *));
    
    for (i = 0; i < Smeta; i++) {
        
        rho[i] = (double *) malloc(Smeta*sizeof(double));
        
        for (j = 0; j < Smeta; j++) {
            
            rho[i][j] = rhovec[k];
            k++;
        }
    }
    
    for (i = 0; i < num_real; i++) {
        
        sample_norep(Slocal, Smeta, local_com_indices);
        
        av_rho[i] = 0.;
        
        for (j = 0; j < Slocal; j++) {
            
            for (k = j+1; k < Slocal; k++) {
                
                av_rho[i] += rho[local_com_indices[j]][local_com_indices[k]];
            }
        }
        
        av_rho[i] *= 2.;
        av_rho[i] /= (double)(Slocal*(Slocal-1));
        
    }
    
    for (i = 0; i < Smeta; i++) {
        
        free(rho[i]);
    }

    free(rho);
    
    free(local_com_indices);
    
    return;
}

void permutations_skew(double *rhovec, int Smeta, int Slocal, int num_real, double *av_rho)
{
    int i, j, k, *local_com_indices = malloc(Slocal*sizeof(int));
    // rho must be the full symmetric metacommunity matrix, expressed as a vector with
    // diagonal elements equal to one
    double **rho, rval;
    time_t tp;
    unsigned long cputime;
    
    time(&tp);
    cputime = (unsigned long)tp;

    srand(cputime);

    k = 0;
    rho = (double **) malloc(Smeta*sizeof(double *));
    
    for (i = 0; i < Smeta; i++) {
        
        rho[i] = (double *) malloc(Smeta*sizeof(double));
        
        for (j = 0; j < Smeta; j++) {
            
            rho[i][j] = rhovec[k];
            k++;
        }
    }
    
    for (i = 0; i < num_real; i++) {
        
        sample_norep(Slocal, Smeta, local_com_indices);
        
        av_rho[i] = 0.;
        
        for (j = 0; j < Slocal; j++) {
            
            for (k = j+1; k < Slocal; k++) {
                
                rval = rho[local_com_indices[j]][local_com_indices[k]];
                if (rval > 0.)
                    av_rho[i] += rval;
                else
                    av_rho[i] += rho[local_com_indices[k]][local_com_indices[j]];
            }
        }
        
        av_rho[i] *= 2.;
        av_rho[i] /= (double)(Slocal*(Slocal-1));
        
    }
    
    for (i = 0; i < Smeta; i++) {
        
        free(rho[i]);
    }
    
    free(rho);
    
    free(local_com_indices);
    
    return;
}

void permutations_skew_ranges(double *rhovec, double *min_env, double *max_env, int Smeta, int Slocal, double env_local, int num_real, double *av_rho)
{
    int i, j, k, size, *local_com_indices = malloc(Slocal*sizeof(int));
    // rho must be the full symmetric metacommunity matrix, expressed as a vector with
    // diagonal elements equal to one
    double **rho, rval;
    time_t tp;
    unsigned long cputime;
    
    time(&tp);
    cputime = (unsigned long)tp;
    
    srand(cputime);
    
    k = 0;
    rho = (double **) malloc(Smeta*sizeof(double *));
    
    for (i = 0; i < Smeta; i++) {
        
        rho[i] = (double *) malloc(Smeta*sizeof(double));
        
        for (j = 0; j < Smeta; j++) {
            
            rho[i][j] = rhovec[k];
            k++;
        }
    }
    
    for (i = 0; i < num_real; i++) {
        
        size = sample_norep_ranges(min_env, max_env, Slocal, Smeta, env_local, local_com_indices);
        
        av_rho[i] = 0.;

        for (j = 0; j < size; j++) {
            
            for (k = j+1; k < size; k++) {
                
                rval = rho[local_com_indices[j]][local_com_indices[k]];
                if (rval > 0.)
                    av_rho[i] += rval;
                else
                    av_rho[i] += rho[local_com_indices[k]][local_com_indices[j]];
            }
        }
        
        av_rho[i] *= 2.;
        av_rho[i] /= (double)(size*(size-1));
    }
    
    for (i = 0; i < Smeta; i++) {
        
        free(rho[i]);
    }
    
    free(rho);
    
    free(local_com_indices);
    
    return;
}

void permutations_skew_sums(double *rhovec, int Smeta, int Slocal, int num_real, double *av_rho)
{
    int i, j, k, *local_com_indices = malloc(Slocal*sizeof(int));
    // rho must be the full symmetric metacommunity matrix, expressed as a vector with
    // diagonal elements equal to one
    double *rho;
    time_t tp;
    unsigned long cputime;
    
    time(&tp);
    cputime = (unsigned long)tp;
    
    srand(cputime);
    
    k = 0;
    rho = (double *) malloc(Smeta*sizeof(double));
    
    for (i = 0; i < Smeta; i++) {
        
        rho[i] = 0.;
        
        for (j = 0; j < Smeta; j++) {
            
            rho[i] += rhovec[k];
            k++;
        }
        
        rho[i] /= (double)Smeta;
    }
    
    for (i = 0; i < num_real; i++) {
        
        sample_norep(Slocal, Smeta, local_com_indices);
        
        av_rho[i] = 0.;
        
        for (j = 0; j < Slocal; j++) {
            
            av_rho[i] += rho[local_com_indices[j]];
        }
        
        av_rho[i] /= (double)Slocal;
    }
    
    free(rho);
    
    free(local_com_indices);
    
    return;
}

void permutations_skew_ranges_weighted(double *rhovec, double *min_env, double *max_env, int Smeta, int Slocal, double env_local, int num_real, double *av_rho, double *cumfreq)
{
    int i, j, k, size, *local_com_indices = malloc(Slocal*sizeof(int));
    // rho must be the full symmetric metacommunity matrix, expressed as a vector with
    // diagonal elements equal to one
    double **rho, rval;
    time_t tp;
    unsigned long cputime;
    
    time(&tp);
    cputime = (unsigned long)tp;
    
    srand(cputime);
    
    k = 0;
    rho = (double **) malloc(Smeta*sizeof(double *));
    
    for (i = 0; i < Smeta; i++) {
        
        rho[i] = (double *) malloc(Smeta*sizeof(double));
        
        for (j = 0; j < Smeta; j++) {
            
            rho[i][j] = rhovec[k];
            k++;
        }
    }
    
    for (i = 0; i < num_real; i++) {
        
        size = sample_norep_ranges_weighted(min_env, max_env, Slocal, Smeta, env_local, local_com_indices, cumfreq);
        
        av_rho[i] = 0.;
        
        for (j = 0; j < size; j++) {
            
            for (k = j+1; k < size; k++) {
                
                rval = rho[local_com_indices[j]][local_com_indices[k]];
                if (rval > 0.)
                    av_rho[i] += rval;
                else
                    av_rho[i] += rho[local_com_indices[k]][local_com_indices[j]];
            }
        }
        
        av_rho[i] *= 2.;
        av_rho[i] /= (double)(size*(size-1));
    }
    
    for (i = 0; i < Smeta; i++) {
        
        free(rho[i]);
    }
    
    free(rho);
    
    free(local_com_indices);
    
    return;
}

void permutations_skew_weighted(double *rhovec, int Smeta, int Slocal, int num_real, double *av_rho,
                                double *cumfreq)
{
    int i, j, k, *local_com_indices = malloc(Slocal*sizeof(int));
    // rho must be the full symmetric metacommunity matrix, expressed as a vector with
    // diagonal elements equal to one
    double **rho, rval;
    time_t tp;
    unsigned long cputime;
    
    time(&tp);
    cputime = (unsigned long)tp;
    
    srand(cputime);
    
    k = 0;
    rho = (double **) malloc(Smeta*sizeof(double *));
    
    for (i = 0; i < Smeta; i++) {
        
        rho[i] = (double *) malloc(Smeta*sizeof(double));
        
        for (j = 0; j < Smeta; j++) {
            
            rho[i][j] = rhovec[k];
            k++;
        }
    }
    
    for (i = 0; i < num_real; i++) {
        
        sample_norep_weighted(Slocal, Smeta, local_com_indices, cumfreq);
        
        av_rho[i] = 0.;
        
        for (j = 0; j < Slocal; j++) {
            
            for (k = j+1; k < Slocal; k++) {
                
                rval = rho[local_com_indices[j]][local_com_indices[k]];
                if (rval > 0.)
                av_rho[i] += rval;
                else
                av_rho[i] += rho[local_com_indices[k]][local_com_indices[j]];
            }
        }
        
        av_rho[i] *= 2.;
        av_rho[i] /= (double)(Slocal*(Slocal-1));
    }
    
    for (i = 0; i < Smeta; i++) {
        
        free(rho[i]);
    }
    
    free(rho);
    
    free(local_com_indices);
    
    return;
}

void sample_norep(int size, int maxi, int *vector)
{
    int num = 0, flag = 0;
    
    do {
        
        vector[num] = (int)(get_random() * (maxi));
        flag = lookfor(vector, vector[num], num);
        if (!flag) num++;
        
    } while (num < size);
    
    return;
}

int sample_norep_ranges(double *min_env, double *max_env, int size, int maxi, double env_local, int *vector)
{
    int num = 0, flag = 0, iter = 0;
    
    do {
        
        vector[num] = (int)(get_random() * (maxi));

        if (min_env[vector[num]] <= env_local && max_env[vector[num]] >= env_local) {
        
            flag = lookfor(vector, vector[num], num);
            if (!flag) num++;
        }
        
        iter++;
        
    } while (num < size && iter < 2*maxi);
    
    return(num);
}

int sample_norep_ranges_weighted(double *min_env, double *max_env, int size, int maxi, double env_local, int *vector, double *cumfreq)
{
    int num = 0, flag = 0, val, iter = 0;
    double u;
    
    do {
        
        u = get_random();
        val = sample_weighted(u, maxi, cumfreq);
        
        if (min_env[val] <= env_local && max_env[val] >= env_local) {
            
            flag = lookfor(vector, val, num);
            if (!flag) {
                
                vector[num] = val;
                num++;
            }
        }
        
        iter++;
        
    } while (num < size && iter < 2*maxi);
    
    return(num);
}

void sample_norep_weighted(int size, int maxi, int *vector, double *cumfreq)
{
    int num = 0, flag = 0, val;
    double u;
    
    do {
        
        u = get_random();
        val = sample_weighted(u, maxi, cumfreq);
        flag = lookfor(vector, val, num);
        if (!flag) {
            
            vector[num] = val;
            num++;
        }
        
    } while (num < size);
    
    return;
}

int lookfor(int *vector, int value, int size)
{
    int i, flag = 0;
    
    for (i = 0; i < size; i++) {
        
        if (vector[i] == value) {
            
            flag = 1;
            break;
        }
    }
    
    return(flag);
}

int sample_weighted(double unif, int size, double *cumfreq)
{
    int i, val = 0;
    
    for (i = 0; i < size; i++) {
        
        if (cumfreq[i] > unif) {
            
            val = i;
            break;
        }
    }
    
    return(val);
}
