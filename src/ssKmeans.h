#ifndef SSKMEANS_H
#define SSKMEANS_H

int ssKmeans(double **x, 
             int p,
             int n,
             int K,
             int *poseq,
             int B,
             int **negeq,
             int negdim,
             int runs,
             int max_iter,
             double **mu,
             int *id,
             double *misc_double);

void initialization_mu(double **x,
                       int n,
                       int p,
                       int K,
                       int B,
                       int *poseq,
                       double **mu);

void membership(double **x,
                int K,
                int n,
                int p,
                int B,
                int negdim,
                int *poseq,
                double **mu,
                double sigma2,
                int *Z,
                int *head,
                int *to,
                int *next,
                int ncomp,
                int *comp_start,
                int *comp_all,
                int *comp_sort_all,
                int *current_Z);

int recompute_means(double **mu,
                    double **x,
                    int *Z,
                    int *poseq,
                    int K,
                    int n,
                    int p,
                    int B);

void array1to2(int a, int b, double *y, double **x);
void array2to1(int a, int b, double *y, double **x);
void array1to2i(int a, int b, int *y, int **x);
void array2to1i(int a, int b, int *y, int **x);

#endif /* SSKMEANS_H */
