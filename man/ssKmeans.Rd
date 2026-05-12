\name{ssKmeans}
\alias{ssKmeans}

\title{Semi-Supervised K-Means Algorithm with Membership Constraints}

\description{
Runs a semi-supervised version of the K-Means algorithm with optional positive and negative membership constraints. Positive constraints require specified observations to share the same class label, while negative constraints require specified blocks of observations to receive different class labels.
}

\usage{
ssKmeans(X, K, tol = 1e-10, pos.eq = NULL, neg.eq = NULL, r = 3, max.iter = 100)
}

\arguments{
  \item{X}{A numeric data matrix with \eqn{n} observations and \eqn{p} variables.}
  \item{K}{The number of clusters.}
  \item{tol}{The convergence tolerance level.}
  \item{pos.eq}{An optional list of positive constraints. Each element of the list contains the row indices of observations that must be assigned the same class label. If \code{NULL}, each observation is treated as its own singleton block.}
  \item{neg.eq}{An optional symmetric adjacency matrix encoding negative constraints between blocks. A value of 1 indicates that the corresponding pair of blocks must be assigned different class labels, and 0 indicates no negative constraint. If \code{NULL}, no negative constraints are imposed.}
  \item{r}{The number of random restarts.}
  \item{max.iter}{The maximum number of iterations allowed for each restart.}
}

\details{
The function implements a semi-supervised extension of the K-Means algorithm that accommodates hard positive and negative constraints.

Positive constraints are handled by forming blocks of observations that must share a common class label. If no positive constraints are supplied, each observation is treated as a singleton block.

Negative constraints are imposed between blocks through the matrix \code{neg.eq}. These constraints require certain blocks to receive different class labels. The algorithm searches over feasible class assignments that satisfy the specified constraints while minimizing the within-cluster sum of squares objective.

Several random restarts are recommended because K-Means-type algorithms are sensitive to initialization. The default number of restarts is \code{r = 3}.
}

\value{
A list containing:
  \item{id}{The final membership/class label vector of length \eqn{n}.}
  \item{size}{The final cluster sizes.}
  \item{centers}{A matrix of final cluster centers with dimension \eqn{K \times p}.}
  \item{var}{The estimated common variance parameter.}
  \item{tss}{The total sum of squares.}
  \item{wcss}{The within-cluster sum of squares.}
  \item{bcss}{The between-cluster sum of squares.}
}

\references{
Melnykov, I. and Melnykov, V. (2020). A note on the formal implementation of the K-means algorithm with hard positive and negative constraints. \emph{Journal of Classification}, 37, 789--809.

Melnykov, V., Melnykov, I. and Michael, S. (2016). Semi-supervised model-based clustering with positive and negative constraints. \emph{Advances in Data Analysis and Classification}, 10, 327--349.
}

\author{
Emma Johnson
}

\seealso{
\code{\link{Kmeans}}
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


\keyword{cluster}
\keyword{classification}
\keyword{multivariate}