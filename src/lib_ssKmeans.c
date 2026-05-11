
#include "array.h"
#include "ssKmeans.h"

void run_ssKmeans(double *x1, int *poseq, int *negeq1, double *mu1, int *id, int *misc_int, double *misc_double){

	int p, n, K, B;
	int negdim, runs, max_iter;
	int **negeq;
		
	double **x, **mu;
	
	p = misc_int[0];
	n = misc_int[1];
	K = misc_int[2];
	B = misc_int[3];
	negdim = misc_int[4];
	runs = misc_int[5];
	max_iter = misc_int[6];
	
	MAKE_MATRIX(x, n, p);
	MAKE_MATRIX(mu, K, p);
	MAKE_MATRIX(negeq, negdim, 2);
	
	array1to2(n, p, x1, x);
	array1to2(K, p, mu1, mu);
	array1to2i(negdim, 2, negeq1, negeq);
	
	ssKmeans(x, p, n, K, poseq, B, negeq, negdim, runs, max_iter, mu, id, misc_double);

	array2to1(K, p, mu1, mu);
		
	FREE_MATRIX(x);
	FREE_MATRIX(mu);
	FREE_MATRIX(negeq);
		
}







