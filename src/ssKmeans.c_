
#include "array.h"
#include "ssKmeans.h"

#ifdef __HAVE_R_
#include <R.h>
#include <Rmath.h>
#endif

#include <math.h>


// initializes first centoids `mu` (updates input mu[k][j])
void initialization_mu(double **x, 
                       int n, 
                       int p, 
                       int K, 
                       int B, 
                       int *poseq, 
                       double **mu) {
  
  double **block_means, *block_prob;
  int *integers, *sample, *block_size;
  
  double rand_num, cum_prob, rem_prob;
  int b, i, j, k;
  int largest;
  
  MAKE_MATRIX(block_means, B, p); // mu_b
  MAKE_VECTOR(block_size, B); // n_b
  MAKE_VECTOR(block_prob, B); // probability mass of block b
  MAKE_VECTOR(integers, B); // block labels 1..B
  MAKE_VECTOR(sample, K); // k values sampled from 1..B
  
  // initializing block_sizes and block_means
  for (b=0; b<B; b++) {
    block_size[b] = 0;
    for (j=0; j<p; j++) {
      block_means[b][j] = 0.0;
    }
  }
  
  // sums of each obs for each block in poseq
  for (i=0; i<n; i++) {
    int b = poseq[i] - 1;
    block_size[b]++;
    for (j=0; j<p; j++) {
      block_means[b][j] += x[i][j];
    }
  }
  
  // final calculation of the block_means
  for (b=0; b<B; b++) {
    for (j=0; j<p; j++) {
      block_means[b][j] = block_means[b][j] / block_size[b];
    }
  }
  
  // initializing an array to generate a sample without replacement
  for (b=0; b<B; b++) {
    integers[b] = b + 1; // block ids
    block_prob[b] = block_size[b] / (double) n;
  }
  
  // picking the numbers without replacement
  largest = B;
  
  // samples K blocks (not just observations)
  for (k=0; k<K; k++) { 
    
#ifdef __HAVE_R_
    rand_num = runif(0.0, 1.0);
#else
    rand_num = rand() / (double) RAND_MAX;
#endif	
    
    b = 0;
    cum_prob = block_prob[0];
    while (rand_num > cum_prob) {
      b++;
      cum_prob += block_prob[b];
    }
    
    sample[k] = b + 1;
    integers[b] = integers[largest - 1];
    rem_prob = block_prob[b];
    block_prob[b] = block_prob[largest - 1];
    
    largest = largest - 1;
    
    // probabilities sum to 1
    for (b=0; b<largest; b++) {
      block_prob[b] = block_prob[b] / (1.0 - rem_prob);
    }
    
  }
  
  // initialize centroids (mu) using mean vectors of the sampled blocks
  for (k=0; k<K; k++) {
    int b = sample[k] - 1;
    for (j=0; j<p; j++) {
      mu[k][j] = block_means[b][j];
    }
  }
  
  FREE_VECTOR(integers);
  FREE_VECTOR(sample);
  FREE_VECTOR(block_prob);
  FREE_MATRIX(block_means);
  FREE_VECTOR(block_size);
  
}


// Build adjacency list
static void adj_list(int B,
                          int negdim,
                          int **negeq,
                          int *head,
                          int *to,
                          int *next) {
  
  //////////////////////////////////////////////////////////////////////////////
  // Adjacency List for negative constraints
  //////////////////////////////////////////////////////////////////////////////
  
  // both (b_u, b_v) & (b_v, b_u)
  // i.e. 1 undirected edge = 2 directed edges (traversal from either direction)
  
  int ep = 0;
  
  for (int b = 0; b < B; b++) {
    head[b] = -1;
  }
  
  for (int i = 0; i < negdim; i++) {
    int b_u = negeq[i][0] - 1;
    int b_v = negeq[i][1] - 1;
    
    // b_u -> b_v (blocks directed one way)
    to[ep] = b_v;
    next[ep] = head[b_u];
    head[b_u] = ep;
    ep++;
    
    // b_v -> b_u (blocks directed the other way)
    to[ep] = b_u;
    next[ep] = head[b_v];
    head[b_v] = ep;
    ep++;
  }

}


// Gather connected components
static int find_comps(int B,
                            int *head,
                            int *to,
                            int *next,
                            int *visited,
                            int *stack,
                            int *comp,
                            int *comp_sort,
                            int *comp_start,
                            int *comp_all,
                            int *comp_sort_all) {
  
  
  int cc = 0;        // component counter
  int offset = 0;    // position in comp_all / comp_sort_all
  
  // clear visited
  for (int b = 0; b < B; b++) {
    visited[b] = 0;
  }
  
  for (int b = 0; b < B; b++) {
    
    if (visited[b]) {
      continue; // skip blocks already traversed
    }
    
    // use DFS to find nodes in block b's connected component
    int sp = 0; // stack-pointer
    int comp_size = 0; // number of blocks in current component
    
    ////////////////////////////////////////////////////////////////////////////
    // Collecting all blocks in b's connected component
    ////////////////////////////////////////////////////////////////////////////
    
    stack[sp] = b; // pop b into stack
    sp++;
    visited[b] = 1; // mark b as visited
    
    while (sp > 0) {
      sp--;
      int b_u = stack[sp]; // 'popped' block/node
      
      comp[comp_size] = b_u;
      comp_size++;
      
      for (int e = head[b_u]; e != -1; e = next[e]) {
        int b_v = to[e];
        if (visited[b_v] == 0) {
          visited[b_v] = 1;
          stack[sp] = b_v;
          sp++;
        }
      }
    } // now the blocks for this component are identified in `comp`
    
    
    
    
    
    
    // record component start
    comp_start[cc] = offset;
    
    // copy comp[] into flat comp_all[]
    for (int c = 0; c < comp_size; c++) {
      comp_all[offset + c] = comp[c];
    }
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    // Sorting blocks in this component by number of neighbors (descending)
    ////////////////////////////////////////////////////////////////////////////
    
    // use `comp_sort` as a selection sort algorithm imposed on blocks in `comp`
    // descending sort by `max_nbrs`
    
    for (int c = 0; c < comp_size; c++) {
      comp_sort[c] = comp[c];
    }
    
    for (int c = 0; c < comp_size; c++) {
      int b_max = c;
      int max_nbrs = 0; // seeks block with most negative constraints
      int b_c = comp_sort[c];
      
      for (int e = head[b_c]; e != -1; e = next[e]) {
        max_nbrs++;
      }
      
      for (int d = c + 1; d < comp_size; d++) {
        int num_nbrs = 0;
        int b_d = comp_sort[d];
        
        for (int e = head[b_d]; e != -1; e = next[e]) {
          num_nbrs++;
        }
        
        if (num_nbrs > max_nbrs) {
          max_nbrs = num_nbrs;
          b_max = d;
        }
      }
      
      if (b_max != c) {
        int swap = comp_sort[c];
        comp_sort[c] = comp_sort[b_max];
        comp_sort[b_max] = swap;
      }
    }
    
    
    
    
    
    
    // copy sorted order into flat comp_sort_all[]
    for (int c = 0; c < comp_size; c++) {
      comp_sort_all[offset + c] = comp_sort[c];
    }
    
    offset += comp_size;
    cc++;
  }
  
  comp_start[cc] = offset;  // sentinel end
  return cc;
}


// updates Z[b] using positive and negative constraints (log negative)
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
                int *current_Z) {

  for (int b = 0; b < B; b++) {
    current_Z[b] = -1;
  }
  
  
  //////////////////////////////////////////////////////////////////////////////
  // computing squared Euclidean Distances per block (dist_sum)
  //////////////////////////////////////////////////////////////////////////////
  
  double **dist_sum;
  MAKE_MATRIX(dist_sum, B, K);
  
  for (int b = 0; b < B; b++) {
    for (int k = 0; k < K; k++) {
      dist_sum[b][k] = 0.0;
    }
  }
  
  for (int i = 0; i < n; i++) {
    int b = poseq[i] - 1; // R to C
    for (int k = 0; k < K; k++) {
      double dist_sq = 0.0;
      for (int j = 0; j < p; j++) {
        double d = x[i][j] - mu[k][j];
        dist_sq += d * d;
      }
      dist_sum[b][k] += dist_sq;
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////
  // pos-only mode: arg_min_k over dist_sum per block
  //////////////////////////////////////////////////////////////////////////////
  
  if (negdim <= 0) { // i.e. no negative contraints
    for (int b = 0; b < B; b++) {
      int arg_min_k = 0;
      double min_dist = dist_sum[b][0];
      for (int k = 1; k < K; k++) {
        if (dist_sum[b][k] < min_dist) {
          min_dist = dist_sum[b][k];
          arg_min_k = k;
        }
      }
      Z[b] = arg_min_k;
    }
    FREE_MATRIX(dist_sum);
    return;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  // neg mode: arg_max_k over feasible block sums (using DFS / backtracking)
  //////////////////////////////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////////////////////////////
  // For each connected component, compute log_sum then
  // choose a feasible labeling maximizing sum log_sum.
  //////////////////////////////////////////////////////////////////////////////
  
  for (int cc = 0; cc < ncomp; cc++) {
    
    int start = comp_start[cc];
    int end = comp_start[cc + 1];
    int comp_size = end - start;
    
    int *comp, *comp_sort;
    MAKE_VECTOR(comp, comp_size);
    MAKE_VECTOR(comp_sort, comp_size);
    
    for (int c = 0; c < comp_size; c++) {
      comp[c] = comp_all[start + c];
      comp_sort[c] = comp_sort_all[start + c];
    }
  
    
    // for singleton block components
    if (comp_size == 1) {
      int b = comp[0];
      int arg_min_k = 0;
      double min_dist = dist_sum[b][0];
      
      for (int k = 1; k < K; k++) {
        if (dist_sum[b][k] < min_dist) {
          min_dist = dist_sum[b][k];
          arg_min_k = k;
        }
      }
      
      Z[b] = arg_min_k;

      continue;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Finding the Log of the Z[b] objective for feasible labelings (log_sum)
    ////////////////////////////////////////////////////////////////////////////
    
    double **log_sum;
    MAKE_MATRIX(log_sum, comp_size, K); // indexed by component-position c, cluster k
    
    // tentative Z assignments
    for (int c = 0; c < comp_size; c++) {
      int b = comp[c];
      current_Z[b] = -1;
    }
    
    for (int c = 0; c < comp_size; c++) {
      for (int k = 0; k < K; k++) {
        log_sum[c][k] = -INFINITY;
      }
    }
    
    int *next_k;
    double *cum_SSE;
    MAKE_VECTOR(next_k, comp_size + 1);
    MAKE_VECTOR(cum_SSE, comp_size + 1);
    
    for (int c = 0; c <= comp_size; c++) {
      next_k[c] = 0;
      cum_SSE[c] = 0.0;
    }
    
    int step = 0; // depth of DFS
    
    while (step >= 0) {
      
      // step 2: use feasible Z's to begin filling log_sum
      // step 4: after backtracking, re-evaluate log_sum w/ new logw 
      if (step == comp_size) {
        
        double SSE  = cum_SSE[comp_size];
        double logw = -SSE / (2.0 * sigma2);
        
        for (int c = 0; c < comp_size; c++) {
          int b = comp[c];
          int k = current_Z[b]; 
          
          // computes the log of the Z objective (keeping things finite)
          if (!isfinite(log_sum[c][k])) {
            log_sum[c][k] = logw;
          }
          else {
            double largest;
            if (log_sum[c][k] > logw) {
              largest = log_sum[c][k];
            } else {
              largest = logw;
            }
            log_sum[c][k] = largest + log(exp(log_sum[c][k] - largest) + exp(logw - largest));
          }
        }
        
        // step 3: backtrack to find BEST k (not just a feasible one)
        step--;
        if (step >= 0) {
          int b = comp_sort[step];
          current_Z[b] = -1;
        }
        continue; // restart while (step >= 0)...
      }
      
      // step 1: create feasible Z assignments (not necessarily best)
      int b_u = comp_sort[step]; // current focus block
      int tried = 0; // 0 = every k was infeasible or exhausted; 1 = some feasible k found
      
      while (next_k[step] < K) {
        int k = next_k[step];
        next_k[step]++;
        
        int ok = 1;
        
        // check whether negative neighbors are already assigned to the same k (b_u, b_v)
        for (int e = head[b_u]; e != -1; e = next[e]) {
          int b_v = to[e];
          if (current_Z[b_v] == k) {
            ok = 0;
            break;
          }
        }
        
        if (ok == 0) { // if negative neighbor is already assigned k,
          continue; // restart while (next_k ...)
        }
        
        current_Z[b_u] = k;
        cum_SSE[step + 1] = cum_SSE[step] + dist_sum[b_u][k];
        
        next_k[step + 1] = 0;
        step++;
        tried = 1;
        break;
      }
      
      if (tried == 0) {
        current_Z[b_u] = -1;
        next_k[step] = 0;
        step--;
      }
    }
    
    FREE_VECTOR(cum_SSE);
    
    // if log_sum is all -INF for this component, there were no feasible labelings
    // then, keep Z = 0 for blocks in this component
    int any_finite = 0;
    for (int c = 0; c < comp_size; c++) {
      for (int k = 0; k < K; k++) {
        if (isfinite(log_sum[c][k])) { 
          any_finite = 1; // i.e. this log_sum is finite
          break; 
        }
      }
      if (any_finite) {
        break;
      }
    }
    
    if (!any_finite) { 
      for (int c = 0; c < comp_size; c++) {
        int b = comp[c];
        Z[b] = 0;
      }
      
      FREE_VECTOR(comp);
      FREE_VECTOR(comp_sort);
      FREE_VECTOR(next_k);
      FREE_MATRIX(log_sum);
      continue; // restart for (int cc...)
    }

    ////////////////////////////////////////////////////////////////////////////
    // Assigning Z[b] based on arg_max_k(log_sum) for current component
    ////////////////////////////////////////////////////////////////////////////
    
    double *cum_log_sum;
    int *arg_max_k;
    double max_log_sum = -INFINITY;
    MAKE_VECTOR(cum_log_sum, comp_size + 1);
    MAKE_VECTOR(arg_max_k, comp_size);
    
    for (int c = 0; c <= comp_size; c++) {
      next_k[c] = 0;
      cum_log_sum[c] = 0.0;
    }
    
    // clear current_Z again
    for (int c = 0; c < comp_size; c++) {
      int b = comp_sort[c];
      current_Z[b] = -1;
      arg_max_k[c] = -1;
    }
    
    step = 0;
    
    // similar to previous while (step ...)
    // now we assign Z based on arg_max_k
    while (step >= 0) {
      
      // step 2: use feasible Z's to find lsum (log Z objective)
      // step 4: after backtracking, re-evaluate lsum w/ new cum_log_sum 
      if (step == comp_size) {
        
        double lsum = cum_log_sum[comp_size];
        
        if (lsum > max_log_sum) {
          max_log_sum = lsum;
          for (int c = 0; c < comp_size; c++) {
            int b = comp[c];
            arg_max_k[c] = current_Z[b];
          }
        }
        
        // step 3: backtrack to find BEST k (not just a feasible one)
        step--;
        if (step >= 0) {
          int b = comp_sort[step];
          current_Z[b] = -1;
        }
        continue;
      }
      
      // step 1.1: feasible Z assignments 
      int b_u = comp_sort[step]; // current focus block
      int tried = 0;
      
      while (next_k[step] < K) {
        int k = next_k[step];
        next_k[step]++;
        
        int ok = 1;
        
        // check whether negative neighbors are already assigned to the same k (b_u, b_v)
        for (int e = head[b_u]; e != -1; e = next[e]) {
          int b_v = to[e];
          if (current_Z[b_v] == k) {
            ok = 0;
            break;
          }
        }
        
        if (ok == 0) { // if negative neighbor is already assigned k,
          continue; // restart while (next_k ...)
        }
        
        // step 1.2: revert back to comp instead of comp_sort for b_u's index
        // i.e. log_sum is ordered by comp, not comp_sort: log_sum[c_ind][k]
        int c_ind = -1;
        for (int c = 0; c < comp_size; c++) {
          if (comp[c] == b_u) { 
            c_ind = c; 
            break; 
          }
        }
        
        // step 1.3: same as before, swapping cum_SSE for cum_log_sum
        if (!isfinite(log_sum[c_ind][k])) {
          continue;
        }
        
        current_Z[b_u] = k;
        cum_log_sum[step + 1] = cum_log_sum[step] + log_sum[c_ind][k];
        
        next_k[step + 1] = 0;
        step++;
        tried = 1;
        break;
      }
      
      if (tried == 0) {
        current_Z[b_u] = -1;
        next_k[step] = 0;
        step--;
      }
    }
    
    FREE_VECTOR(next_k);
    FREE_VECTOR(cum_log_sum);
    
    // Write output Z for this component w/ arg_max_k
    for (int c = 0; c < comp_size; c++) {
      int b = comp[c];
      
      if (arg_max_k[c] == -1) {
        Z[b] = 0;
      } else {
        Z[b] = arg_max_k[c];
      }
    }
    
    FREE_VECTOR(comp);
    FREE_VECTOR(comp_sort);
    FREE_VECTOR(arg_max_k);
    FREE_MATRIX(log_sum);
  } // end for (int cc...)
  
  FREE_MATRIX(dist_sum);
} 


// recomputes `mu` based on membership assignments (updates input mu[k][j])
int recompute_means(double **mu, 
                    double **x, 
                    int *Z, 
                    int *poseq, 
                    int K, 
                    int n, 
                    int p, 
                    int B) {
  
  int i, k, j = 0;
  int flag; // i.e. stopper
  int *nk;
  
  MAKE_VECTOR(nk, K); // kth cluster size
  
  flag = 0;
  
  // re-initializing mu full of 0s
  for (k=0; k<K; k++) {
    nk[k] = 0;
    for (j=0; j<p; j++) {
      mu[k][j] = 0.0;
    }
  }
  
  // recompute centroids based on membership assignments
  for (i=0; i<n; i++) {
    int b = poseq[i] - 1;
    int k = Z[b];
    for (j=0; j<p; j++) {
      mu[k][j] += x[i][j];
    }
    nk[k]++;
  }
  
  // centroids heading into next iteration
  for (k=0; k<K; k++) {
    if (nk[k] != 0) {
      for (j=0; j<p; j++) {
        mu[k][j] = mu[k][j] / nk[k];
      }
    } else {
      flag = 1;
      break;
    }
  }
  
  FREE_VECTOR(nk);
  
  return flag;
  // 0 : success (all clusters non-empty)
  // 1 : failure (at least one empty cluster)
  
}


// main shell that calls above 3 functions
int ssKmeans(double **x, 
             int p, 
             int n, 
             int K, 
             int *poseq, 
             int  B, 
             int **negeq, 
             int negdim, 
             int runs, 
             int max_iter, 
             double **mu, 
             int *id, 
             double *misc_double) {
  
  double sigma2, best_sigma2;
  double wss, old_wss, best_wss;
  double **best_mu;
  double dist_sq, min_dist_sq, eps;

  int i, j, k, l;
  int *Z;
  int flag, iter;
  
  eps = misc_double[0]; // i.e. tol
  
  MAKE_MATRIX(best_mu, K, p);
  best_wss = INFINITY;
  best_sigma2 = INFINITY;
  
  MAKE_VECTOR(Z, B);
  
  int *head, *to, *next;
  head = NULL;
  to = NULL;
  next = NULL;
  
  int negdim2;
  
  if (negdim > 0) {
    negdim2 = 2 * negdim;
    MAKE_VECTOR(head, B);
    MAKE_VECTOR(to, negdim2);
    MAKE_VECTOR(next, negdim2);
    adj_list(B, negdim, negeq, head, to, next);
  } else {
    negdim2 = 0;
  }
  
  int *visited, *stack, *comp, *comp_sort, *comp_start, *comp_all, *comp_sort_all;
  
  MAKE_VECTOR(visited, B);
  MAKE_VECTOR(stack, B);
  MAKE_VECTOR(comp, B);
  MAKE_VECTOR(comp_sort, B);
  MAKE_VECTOR(comp_start, B + 1);
  MAKE_VECTOR(comp_all, B);
  MAKE_VECTOR(comp_sort_all, B);
  
  int ncomp = 0;
  
  if (negdim > 0) {
    ncomp = find_comps(B, head, to, next, visited, stack, comp, comp_sort,
                             comp_start, comp_all, comp_sort_all);
  }
  
  int *current_Z;
  MAKE_VECTOR(current_Z, B);
  
  
  for (l=0; l<runs; l++) { // number of Kmeans restarts (runs = nstart)
    
    old_wss = INFINITY;
    iter = 1; // variable for keeping track of # of Kmeans restarts
    
    // Obtain initial cluster centers
    initialization_mu(x, n, p, K, B, poseq, mu);
    
    // Estimate sigma based on distances to current cluster centers
    wss = 0.0;
    for (i=0; i<n; i++) {
      min_dist_sq = INFINITY;
      for (k=0; k<K; k++) {
        dist_sq = 0.0;
        for (j=0; j<p; j++) { 
          double dist = x[i][j] - mu[k][j];
          dist_sq += dist * dist;
        }
        if (dist_sq < min_dist_sq) {
          min_dist_sq = dist_sq;
        }
      }
      wss += min_dist_sq;
    }
    
    // initialize variance 
    sigma2 = wss / ((double)p *(double) n);
    
    flag = 0;
    
    // fabs = floating point absolute value
    while ((flag == 0) && (iter < max_iter) && (fabs((wss - old_wss) / wss) > eps)) {
      
      membership(x, K, n, p, B, negdim, poseq, mu, sigma2, Z, head, to, next,
                 ncomp, comp_start, comp_all, comp_sort_all, current_Z);
      
      
      flag = recompute_means(mu, x, Z, poseq, K, n, p, B); // flag = 1 means cluster is completely dissolved
      
      old_wss = wss;
      
      wss = 0.0;
      for (i=0; i<n; i++) {
        int b = poseq[i] - 1;
        int k = Z[b];
        for (j=0; j<p; j++) {
          double dist = x[i][j] - mu[k][j];
          wss += dist * dist;
        }
      }
      
      sigma2 = wss / ((double)p * ((double)n - (double)K));
      
      iter++;
      
    } // end of while
    
    if (wss < best_wss) { 
      best_wss = wss;
      best_sigma2 = sigma2;
      for (i=0; i<n; i++) {
        int b = poseq[i] - 1;
        id[i] = Z[b];
      }
      for (k=0; k<K; k++) {
        for (j=0; j<p; j++) {
          best_mu[k][j] = mu[k][j];
        }
      }
    }
    
  } // end of loop over l
  
  for (k=0; k<K; k++) {
    for (j=0; j<p; j++) {
      mu[k][j] = best_mu[k][j];
    }
  }
  misc_double[1] = best_sigma2;
  misc_double[2] = best_wss;
  
  FREE_VECTOR(Z);
  
  if (negdim > 0) {
    FREE_VECTOR(next);
    FREE_VECTOR(to);
    FREE_VECTOR(head);
  }
  
  FREE_MATRIX(best_mu);
  FREE_VECTOR(comp_sort_all);
  FREE_VECTOR(comp_all);
  FREE_VECTOR(comp_start);
  FREE_VECTOR(comp_sort);
  FREE_VECTOR(comp);
  FREE_VECTOR(stack);
  FREE_VECTOR(visited);
  
  FREE_VECTOR(current_Z);
  
  
  return 0;
  
}
