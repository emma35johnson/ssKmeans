# ssKmeans

`ssKmeans` is an R package developed as part of my Master's thesis, *Accommodating Positive & Negative Constraints in the K-Means Algorithm*. The package implements a semi-supervised extension of the K-Means algorithm that accommodates hard positive and negative membership constraints.

Standard K-Means is an unsupervised clustering algorithm that partitions observations into `K` clusters by minimizing the within-cluster sum of squares. In practice, partial prior information may be available. For example, some observations may be known to belong together, while others may be known to belong to different groups. The `ssKmeans` package incorporates this information through positive and negative constraints.

Positive constraints require specified observations to receive the same class label. These constraints are handled through block formation. Negative constraints require specified observations or blocks of observations to receive different class labels. These constraints are handled using graph-based structures and feasible class assignment enumeration. The user-facing interface is written in R, while the computational backend is implemented in C.

## Installation

From the package directory, run:

```r
devtools::install()
```

Then load the package:

```r
library(ssKmeans)
```

## Main Function

```r
ssKmeans(X, K, tol = 1e-10, pos.eq = NULL, neg.eq = NULL,
         neg.type = c("matrix", "edge"), r = 3, max.iter = 100)
```

## Basic Example

```r
library(ssKmeans)

X <- as.matrix(iris[, 1:4])
K <- 3

set.seed(5)

ss.none <- ssKmeans(
  X = X,
  K = K,
  r = 50
)

ss.none
```

## Positive Constraints

Positive constraints are supplied as a list. Each element of the list contains the row indices of observations that must receive the same class label.

```r
actual <- as.numeric(iris$Species)

set.seed(5)

pos <- lapply(
  split(seq(actual), actual),
  sample,
  size = 5
)

ss.pos <- ssKmeans(
  X = X,
  K = K,
  pos.eq = pos,
  r = 50
)

ss.pos
```

## Negative Constraints

Negative constraints can be supplied as an adjacency matrix or as an edge list. Below is the edge list representation, where each row gives a pair of observations that must receive different class labels.

```r
neg.edge <- rbind(
  c(42, 58),
  c(51, 110),
  c(99, 107),
  c(99, 110)
)

ss.neg <- ssKmeans(
  X = X,
  K = K,
  neg.eq = neg.edge,
  neg.type = "edge",
  r = 50
)

ss.neg
```

## Positive and Negative Constraints Together

When positive constraints are supplied, observations are first grouped into blocks. Negative constraints can then be imposed between blocks.

```r
set.seed(5)

pos <- lapply(
  split(seq(actual), actual),
  sample,
  size = 5
)

neg.block <- rbind(
  c(1, 2),
  c(1, 3),
  c(2, 3)
)

ss.both <- ssKmeans(
  X = X,
  K = K,
  pos.eq = pos,
  neg.eq = neg.block,
  neg.type = "edge",
  r = 50
)

ss.both
```

## Returned Values

The function returns a list containing:

- `id`: final membership/class label vector
- `size`: final cluster sizes
- `centers`: final cluster centers
- `var`: estimated common variance parameter
- `tss`: total sum of squares
- `wcss`: within-cluster sum of squares
- `bcss`: between-cluster sum of squares

## Thesis Context

This package accompanies the thesis Accommodating Positive & Negative Constraints in the K-Means Algorithm. The method preserves the iterative structure and within-cluster sum of squares objective of classical K-Means while incorporating prior information through hard membership constraints.
