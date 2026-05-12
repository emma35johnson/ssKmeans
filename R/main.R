
ssKmeans <- function(X, K, tol = 1e-10, pos.eq = NULL, neg.eq = NULL,
                     neg.type = c("matrix", "edge"),
                     r = 3, max.iter = 100){
  
  if (K < 1) stop("Wrong number of mixture components K...\n")
  if (tol <= 0) stop("Wrong value of tol...\n")
  if (r < 1) stop("Wrong number of random restarts r...\n")
  if (max.iter < 1) stop("Wrong number of iterations iter...\n")
  
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  x1 <- as.vector(t(X))
  id <- rep(0, n)
  
  if (is.null(pos.eq)){
    poseq1 <- 1:n
  } else {
      if (!is.list(pos.eq)){
        stop("`pos.eq` must be a block-indexed list of row indices.")
      }
    B_input <- length(pos.eq)
    poseq1 <- integer(n)
    for (b in 1:B_input){
      idx <- as.integer(pos.eq[[b]])
      poseq1[idx] <- b
    }
    remaining <- which(poseq1 == 0L)
    if (length(remaining) > 0L) {
      poseq1[remaining] <- seq(from = B_input + 1L, length.out = length(remaining))
      }
    }
  B <- max(poseq1)
  
  if (is.null(neg.eq)){
    negdim <- 0L
    negeq1 <- integer(0)
  } else {
    if (neg.type == "matrix") {
      neg_mat <- as.matrix(neg.eq)
      if (nrow(neg_mat) != ncol(neg_mat)){
        stop("When `neg.type = 'matrix'`, `neg.eq` must be a square adjacency matrix.")
      }
      diag(neg_mat) <- 0
      edges <- which(neg_mat != 0 & upper.tri(neg_mat), arr.ind = TRUE)
    } else if (neg.type == "list") {
      edges <- as.matrix(neg.eq)
      if (ncol(edges) != 2){
        stop("When `neg.type = 'list'`, `neg.eq` must have two columns.")
      }
      edges <- matrix(as.integer(edges), ncol = 2)
    }
    negdim <- nrow(edges)
    if (negdim == 0L){
      negeq1 <- 0L
    } else {
      negeq1 <- as.vector(t(edges))
    }
  }
  
  mu1 <- rep(0, K*p)
  
  misc_int <- c(p, n, K, B, negdim, r, max.iter)
  misc_double <- c(tol, 0.0, 0.0)
  
  Q <- .C("run_ssKmeans", x1 = as.double(x1), poseq = as.integer(poseq1), negeq1 = as.integer(negeq1), mu1 = as.double(mu1), id = as.integer(id), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), PACKAGE = "ssKmeans")
  
  ssb <- ifelse(K == 1, 0, sum(diag(var(X) * (n-1))) - Q$misc_double[3])
  
  ret <- list(id = Q$id + 1, 
              size = tabulate(Q$id + 1, nbins = K),
              centers = matrix(Q$mu1, nrow = K, byrow = TRUE), 
              var = Q$misc_double[2], 
              tss = Q$misc_double[3] + ssb,
              wcss = Q$misc_double[3], 
              bcss = ssb)
  
  class(ret) <- "ssKmeans"
  return(ret)
  
}

