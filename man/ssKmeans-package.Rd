\name{ssKmeans-package}
\alias{ssKmeans-package}
\alias{ssKmeans-package}
\docType{package}
\title{
Semi-supervised Kmeans algorithm with positive and negative constraints.
}
\description{
The package runs the Kmeans algorithm with positive and negative constraints provided. Positive constraints assign data points to belong to the same cluster. Negative restrictions require observations to belong to different clusters.
}
\details{
\tabular{ll}{
Package: \tab ssKmeans\cr
Type: \tab Package\cr
Version: \tab 0.1.1\cr
Date: \tab 2026-05-19\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab no\cr
}
Function 'ssKmeans' runs the Kmeans algorithm with constraints.
}
\author{ Volodymyr Melnykov and Igor Melnykov

Maintainer: Volodymyr Melnykov <vmelnykov@cba.ua.edu>
}
\references{
Melnykov, V., Melnykov, I. and Michael, S. (2015) Semi-supervised model-based clustering with positive and negative constraints, accepted by Advances in Data Analysis and Classification
}
\examples{

set.seed(123)

X <- as.matrix(iris[,-5])

# K-means with the last 50 points required to belong to the same cluster

peq <- c(1:100, rep(101, 50))
Q <- ssKmeans(X, K = 3, pos.eq = peq); Q$id

}
