
# Iris illustration for ssKmeans 
# (shown in thesis)

library(ssKmeans)
library(mclust)

# Data setup

X <- as.matrix(iris[, 1:4])
actual <- as.numeric(iris$Species)
K <- 3
n <- nrow(X)

################################################################################

# Baseline kmeans

set.seed(5)
km <- kmeans(X, centers = K, nstart = 50)
adjustedRandIndex(km$cluster, actual) # ARI
km$tot.withinss # WCSS
km$betweenss # BCSS
km$tot.withinss / (nrow(X) * ncol(X)) # variance

################################################################################

# helper fxn for ssKmeans results

ssRes <- function(ss.res) {
  list(
    ARI = adjustedRandIndex(ss.res$id, actual) ,
    WCSS = ss.res$wcss,
    BCSS = ss.res$bcss,
    Var = ss.res$var
  )
}

################################################################################

# ssKmeans with no constraints

set.seed(5)
ss.none <- ssKmeans(
  X = X,
  K = K,
  r = 50
)
ssRes(ss.none)

################################################################################

# ssKmeans w/ positive-only constraints

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
ssRes(ss.pos)

################################################################################

# ssKmeans w/ negative-only constraints

set.seed(5)
y <- lapply(pos, sample, size = 2)
neg <- rbind(
  cbind(y[[1]], y[[2]]),
  cbind(y[[1]], y[[3]]),
  cbind(y[[2]], y[[3]])
)

set.seed(5)
ss.neg <- ssKmeans(
  X = X,
  K = K,
  neg.eq = neg,
  neg.type = "list",
  r = 50
)
ssRes(ss.neg)


################################################################################

# ssKmeans w/ positive-negative constraints

neg2 <- rbind(
  c(1, 2),
  c(1, 3),
  c(2, 3)
)

set.seed(5)
ss.both <- ssKmeans(
  X = X,
  K = K,
  pos.eq = pos,
  neg.eq = neg2,
  neg.type = "list",
  r = 50
)
ssRes(ss.both)

################################################################################

# ssKmeans w/ conflicting positive constraints

set.seed(5)
pos.conflict <- list(
  c(sample(1:50, size = 1),
    sample(51:100, size = 1),
    sample(101:150, size = 1))
  )

set.seed(5)
ss.conflict1 <- ssKmeans(
  X = X,
  K = K,
  pos.eq = pos.conflict,
  r = 50
)
ssRes(ss.conflict1)

################################################################################

# ssKmeans w/ conflicting negative constraints

set.seed(5)
neg.conflict <- do.call(rbind, lapply(pos, sample, size = 2))

set.seed(5)
ss.conflict2 <- ssKmeans(
  X = X,
  K = K,
  neg.eq = neg.conflict,
  neg.type = "list",
  r = 50
)
ssRes(ss.conflict2)
