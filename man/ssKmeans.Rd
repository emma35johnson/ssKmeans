\name{ssKmeans}
\alias{ssKmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Semi-supervised Kmeans algorithm with membership constraints }
\description{
Runs the Kmeans algorithm with positive and negative constraints provided. Positive constraints indicate the points that have to belong to the  same cluster. Negative constraints are used to restrict observations to belong to different clusters.
}
\usage{
ssKmeans(X, K, tol = 1e-10, pos.eq = NULL, neg.eq = NULL, r = 3, max.iter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset matrix (n x p) }
  \item{K }{ number of clusters }
  \item{tol }{ tolerance level }
  \item{pos.eq }{ set of positive constraints (vector of length n)}
  \item{neg.eq }{ set of negative constraints (m x 2 matrix of pairwise constraints)}
  \item{r }{ number of restarts}
  \item{max.iter }{ maximum number of iterations allowed}
}
\details{
Runs the Kmeans algorithm with positive and negative equivalence constraints. The function returns estimated cluster centers, partitioning, variance estimate, and within- and between-cluster sums of squares. In the vector 'pos.eq', observations that have to be in the same cluster, must belong to the same block. If positive equivalences are not provided, each observations is treated as a separate block of size 1. In this case, 'pos.eq' is a sequence of integers from 1 to n. Negative constraints must be organized in the form of a m x 2 matrix, where m is the number of pairwise restrictions for all blocks.
To provide better performance of the algorithm, several restarts are recommended. The default number of restarts is 3.
}
\value{
  \item{id }{ classification vector (length n)}
  \item{centers }{ matrix of cluster centers (K x p)}
  \item{var }{ variance}
  \item{ssw }{ within-cluster sum of squares}
  \item{ssb }{ between-cluster sum of squares}
}
\references{
Melnykov, V., Melnykov, I. and Michael, S. (2015) Semi-supervised model-based clustering with positive and negative constraints, accepted by Advances in Data Analysis and Classification}
\author{ Melnykov, V. and Melnykov, I.}

\seealso{kmeans}
\examples{

set.seed(123)

X <- as.matrix(iris[,-5])

# Traditional Kmeans without any constraints

Q1 <- ssKmeans(X, K = 3); Q1$id


# K-means with the last 50 points required to belong to the same cluster

peq <- c(1:100, rep(101, 50))
Q2 <- ssKmeans(X, K = 3, pos.eq = peq); Q2$id


# K-means with the points 131 to 140 in one cluster and points 141:150 in another cluster

peq <- c(1:130, rep(131, 10), rep(132, 10))
neq <- rbind(c(131, 132), c(132, 131))
Q3 <- ssKmeans(X, 3, pos.eq = peq, neg.eq = neq); Q3$id

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Kmeans }
\keyword{ Semi-supervised clustering }
\keyword{ Positive and negative constraints }% __ONLY ONE__ keyword per line
