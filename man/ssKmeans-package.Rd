\name{ssKmeans-package}
\alias{ssKmeans-package}
\docType{package}

\title{
Semi-Supervised K-Means Algorithm with Positive and Negative Constraints
}

\description{
The \pkg{ssKmeans} package implements a semi-supervised version of the K-Means algorithm with optional hard positive and negative membership constraints. Positive constraints require specified observations to share the same class label, while negative constraints require specified blocks of observations to receive different class labels.
}

\details{
\tabular{ll}{
Package: \tab ssKmeans\cr
Type: \tab Package\cr
Version: \tab 0.1.1\cr
Date: \tab 2026-05-12\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab no\cr
}

The main function, \code{\link{ssKmeans}}, runs a semi-supervised K-Means algorithm that accommodates hard positive and negative constraints. Positive constraints are used to form blocks of observations that must receive a common class label. Negative constraints are then imposed between blocks to require different class labels. The method returns the final membership labels, cluster centers, cluster sizes, variance estimate, and within- and between-cluster sums of squares.
}

\author{
Emma Johnson

Maintainer: Emma Johnson
}

\references{
Melnykov, I. and Melnykov, V. (2020). A note on the formal implementation of the K-means algorithm with hard positive and negative constraints. \emph{Journal of Classification}, 37, 789--809.

Melnykov, V., Melnykov, I. and Michael, S. (2016). Semi-supervised model-based clustering with positive and negative constraints. \emph{Advances in Data Analysis and Classification}, 10, 327--349.
}

\seealso{
\code{\link{ssKmeans}},
\code{\link[stats]{kmeans}}
}

\examples{
  set.seed(5)
  
  X <- as.matrix(iris[, 1:4])
  species <- as.numeric(iris$Species)
  K <- 3
  
  # ssKmeans with no constraints
  # (equiv. to standard kmeans)
  
  ss.none <- ssKmeans(
    X = X,
    K = K,
    r = 50
  )
  
  ss.none
  
  
  # ssKmeans with positive-only constraints
  # (blocks of observations)
  
  pos <- list(
    c(16, 42), 
    c(51, 61, 71), 
    c(106, 107, 132)
  )
  
  pos
  
  ss.pos <- ssKmeans(
    X = X,
    K = K,
    pos.eq = pos,
    r = 50
  )
  
  ss.pos
  
  
  # ssKmeans with negative-only constraints 
  # (observation-level)
  
  neg.edge <- rbind(
    c(42, 58),
    c(51, 110),
    c(99, 107),
    c(99, 110)
  )
  
  neg.edge
  
  ss.neg <- ssKmeans(
    X = X,
    K = K,
    neg.eq = neg.edge,
    neg.type = "edge",
    r = 50
  )
  
  ss.neg
  
  
  # ssKmeans with both positive and negative constraints
  # (block-level)
  
  neg.mat <- matrix(
    rep(1, 9),
    nrow = 3
  ) 
  diag(neg.mat) <- rep(0, 3)
  
  neg.mat
  
  ss.both <- ssKmeans(
    X = X,
    K = K,
    pos.eq = pos,
    neg.eq = neg.mat,
    neg.type = "matrix",
    r = 50
  )
  
  ss.both
}


\keyword{package}
\keyword{cluster}
\keyword{classification}
\keyword{multivariate}