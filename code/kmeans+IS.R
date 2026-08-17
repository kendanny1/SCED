### ============================================================
### Packages
### ============================================================
library(expm)       # Matrix exponential and related matrix operations
library(openxlsx)   # Reading and writing Excel files
library(clue)       # Linear assignment problem (solve_LSAP)

### ============================================================
### Global parameters
### ============================================================
p = 15              # Dimension of each observation
KC = 7              # Number of clusters
jj = 6              # Simulation / dataset index
n = 250             # Number of observations in the current dataset


### ============================================================
### Adjusted Rand Index (ARI)
### ============================================================
# Compute the Adjusted Rand Index (ARI) between two clustering
# assignments. ARI measures the agreement between two partitions
# while correcting for agreement that could occur by chance.
#
# x: clustering assignment vector
# y: reference / true clustering assignment vector
#
# The function is equivalent to the standard ARI definition based
# on the contingency table between x and y.
adjustedRandIndex = function(x, y)
{
  x <- as.vector(x)
  y <- as.vector(y)
  
  if (length(x) != length(y))
    stop("arguments must be vectors of the same length")
  
  # Contingency table between the two clustering assignments
  tab <- table(x, y)
  
  # If both partitions contain only one cluster, they are identical
  if (all(dim(tab) == c(1, 1)))
    return(1)
  
  # Number of pairs assigned to the same cluster in both partitions
  a <- sum(choose(tab, 2))
  
  # Number of pairs assigned to the same cluster in x
  # but to different clusters in y
  b <- sum(choose(rowSums(tab), 2)) - a
  
  # Number of pairs assigned to different clusters in x
  # but to the same cluster in y
  c <- sum(choose(colSums(tab), 2)) - a
  
  # Number of pairs assigned to different clusters
  # in both partitions
  d <- choose(sum(tab), 2) - a - b - c
  
  # Adjusted Rand Index
  ARI <- (
    a - (a + b) * (a + c) / (a + b + c + d)
  ) / (
    (a + b + a + c) / 2 -
      (a + b) * (a + c) / (a + b + c + d)
  )
  
  return(ARI)
}


### ============================================================
### Load data
### ============================================================

# Load the simulated dataset corresponding to simulation index jj.
# The file is expected to contain:
#   Xlg : observations
#   Clg : true cluster labels
#   kk  : current dataset / replicate index
#   beta0, mu1, ..., mu7 : true parameter values
load(
  paste(
    "/Users/kendanny/Desktop/R scripts/sep/kotzp15/data_",
    jj,
    ".RData",
    sep = ""
  )
)

# Extract the n observations corresponding to the current replicate.
X = Xlg[((kk - 1) * n + 1):(kk * n),]


### ============================================================
### Standard K-means clustering
### ============================================================

# Apply K-means clustering using Lloyd's algorithm.
#
# iter.max specifies the maximum number of iterations allowed
# for the K-means algorithm.
kx = kmeans(
  X,
  centers = KC,
  algorithm = "Lloyd",
  iter.max = 100
)

# Estimated cluster labels from K-means
cx = kx$cluster

# Store the clustering result
Cle = c(Cle, cx)


### ============================================================
### Clustering performance: Adjusted Rand Index
### ============================================================

# Compare the estimated clustering labels with the true labels.
ri = adjustedRandIndex(
  cx,
  Clg[((kk - 1) * n + 1):(kk * n)]
)

# Store ARI for the current replicate
rii = c(rii, ri)


### ============================================================
### Estimate cluster means
### ============================================================

# mu is a p x KC matrix.
# Column i contains the estimated mean vector of cluster i.
mu = matrix(0, p, KC)

for (i in 1:KC) {
  
  # If the cluster contains only one observation, directly use
  # that observation as the cluster mean.
  #
  # This special case is needed because X[cx == i, ] may be
  # returned as a vector rather than a matrix when the cluster
  # contains only one observation.
  if (length(X[cx == i, ]) == p) {
    
    mu[, i] = X[cx == i, ]
    
  } else {
    
    # Compute the sample mean of cluster i
    mu[, i] =
      colSums(X[cx == i, ]) / dim(X[cx == i, ])[1]
  }
}


### ============================================================
### Mean estimation error
### ============================================================

# The cluster labels produced by K-means are arbitrary.
# For example, estimated cluster 1 does not necessarily correspond
# to true cluster 1. Therefore, we first solve a linear assignment
# problem to find the optimal matching between estimated and true
# cluster means.

# True cluster means
mu_true = cbind(
  mu1, mu2, mu3, mu4, mu5, mu6, mu7
)

# Estimated cluster means
mu_est = mu[, 1:7]

# cost[i, j] = squared Euclidean distance between
# estimated cluster mean i and true cluster mean j.
cost = matrix(0, KC, KC)

for (i in 1:KC) {
  for (j in 1:KC) {
    
    cost[i, j] =
      sum((mu_est[, i] - mu_true[, j])^2)
  }
}

# Find the one-to-one assignment that minimizes the
# total squared distance between estimated and true means.
assignment = solve_LSAP(cost)

# Minimum total squared error after optimal label matching
min_sse = sum(
  cost[cbind(1:KC, assignment)]
)

# Root mean squared error of the estimated cluster means.
#
# Divide by KC * p because there are KC mean vectors,
# each containing p parameters.
MSEgwkm = c(
  MSEgwkm,
  sqrt(min_sse / (KC * p))
)


### ============================================================
### Estimate the common covariance matrix
### ============================================================

# True covariance matrix used in the simulation.
#
# The covariance structure has:
#   - off-diagonal elements = 0.3
#   - diagonal elements     = 1.0
#
# Therefore:
# Sigma = 0.3 * 11' + 0.7 * I_p
#
# Note that the expression below generates a matrix with
# diagonal = 0.3 + 0.7 = 1 and off-diagonal = 0.3.
sig = matrix(0.3, p, p) + 0.7 * diag(p)

# Initialize the estimated covariance matrix
Csig = matrix(0, p, p)

for (i in 1:KC) {
  
  # Calculate the within-cluster scatter matrix for cluster i.
  #
  # Again, handle the one-observation case separately because
  # R drops dimensions when subsetting a matrix to one row.
  if (length(X[cx == i, ]) == p) {
    
    # Outer product of the centered observation
    Csig =
      Csig +
      outer(
        X[cx == i, ] - mu[, i],
        X[cx == i, ] - mu[, i],
        FUN = "*"
      )
    
  } else {
    
    # Cross-product of the centered observations
    Csig =
      Csig +
      (t(X[cx == i, ]) - mu[, i]) %*%
      t(t(X[cx == i, ]) - mu[, i])
  }
}

# Convert the total within-cluster scatter matrix into
# the empirical covariance estimate.
Csig = Csig / n


### ============================================================
### Covariance estimation error
### ============================================================

# Compare the estimated covariance matrix with the true covariance
# matrix. Since the covariance matrix is symmetric, only the
# upper-triangular elements (including the diagonal) are used.
#
# There are p(p+1)/2 unique elements in a symmetric p x p matrix.
MSEkmsig = c(
  MSEkmsig,
  sqrt(
    sum(
      (
        Csig[upper.tri(Csig, diag = TRUE)] -
          jj^2 / 100 *
          sig[upper.tri(sig, diag = TRUE)]
      )^2
    ) /
      (p * (p + 1) / 2)
  )
)


### ============================================================
### Iterative reassignment using the estimated covariance matrix
### ============================================================

# Initialize the cluster assignments.
# cx2 contains the updated cluster labels.
cx2 = rep(0, n)

# Continue updating cluster assignments until they converge.
#
# The algorithm iteratively performs:
#   1. Compute the Mahalanobis distance from each observation
#      to each cluster mean.
#   2. Assign each observation to the nearest cluster.
#   3. Re-estimate cluster means.
#   4. Re-estimate the common covariance matrix.
#
# The loop terminates when the cluster assignments no longer change.
while (sum(cx == cx2) < n) {
  
  # Save the previous cluster assignments.
  cx = cx2
  
  # Reset the updated cluster assignment vector.
  cx2 = c()
  
  # Yu[i, j] stores the squared Mahalanobis distance between
  # observation i and cluster j.
  Yu = matrix(0, n, KC)
  
  for (i in 1:KC) {
    
    # Squared Mahalanobis distance:
    #
    #   (x - mu_i)' Sigma^{-1} (x - mu_i)
    #
    # Here, the calculation is vectorized over all n observations.
    Yu[, i] =
      colSums(
        (
          solve(Csig) %*%
            (t(X) - mu[, i])
        ) *
          (t(X) - mu[, i])
      )
  }
  
  # Assign each observation to the cluster with the smallest
  # squared Mahalanobis distance.
  for (i in 1:n) {
    cx2 = c(
      cx2,
      which.min(Yu[i, ])
    )
  }
  
  
  ### ----------------------------------------------------------
  ### Update cluster means
  ### ----------------------------------------------------------
  
  for (i in 1:KC) {
    
    # Handle the case where cluster i contains only one observation.
    if (length(X[cx2 == i, ]) == p) {
      
      mu[, i] = X[cx2 == i, ]
      
    } else {
      
      # Re-estimate the mean of cluster i.
      mu[, i] =
        colSums(X[cx2 == i, ]) /
        dim(X[cx2 == i, ])[1]
    }
  }
  
  
  ### ----------------------------------------------------------
  ### Update covariance matrix
  ### ----------------------------------------------------------
  
  Csig = matrix(0, p, p)
  
  for (i in 1:KC) {
    
    if (length(X[cx2 == i, ]) == p) {
      
      Csig =
        Csig +
        outer(
          X[cx2 == i, ] - mu[, i],
          X[cx2 == i, ] - mu[, i],
          FUN = "*"
        )
      
    } else {
      
      Csig =
        Csig +
        (t(X[cx2 == i, ]) - mu[, i]) %*%
        t(t(X[cx2 == i, ]) - mu[, i])
    }
  }
  
  # Estimate the common covariance matrix using all observations.
  Csig = Csig / n
}


### ============================================================
### Store final clustering result
### ============================================================

# Store the final cluster assignment obtained from
# the iterative Mahalanobis-distance procedure.
Clek = c(Clek, cx2)


### ============================================================
### Clustering performance after iterative reassignment
### ============================================================

# Compute ARI using the final cluster assignment.
ri2 = adjustedRandIndex(
  cx,
  Clg[((kk - 1) * n + 1):(kk * n)]
)


### ============================================================
### Observation-level mean estimation error
### ============================================================

# Construct the estimated mean vector corresponding to each
# observation according to its final cluster assignment.
betahe = c()

for (i in 1:n) {
  betahe = c(
    betahe,
    mu[, cx2[i]]
  )
}

# Compare the estimated cluster mean assigned to each observation
# with its true mean vector beta0.
MSEhe = c(
  MSEhe,
  sqrt(
    sum((betahe - beta0)^2) / (n * p)
  )
)


### ============================================================
### Cluster mean estimation error after iterative reassignment
### ============================================================

# True cluster means
mu_true = cbind(
  mu1, mu2, mu3, mu4, mu5, mu6, mu7
)

# Estimated cluster means
mu_est = mu[, 1:7]

# Construct the cost matrix for optimal label matching.
cost = matrix(0, KC, KC)

for (i in 1:KC) {
  for (j in 1:KC) {
    
    cost[i, j] =
      sum((mu_est[, i] - mu_true[, j])^2)
  }
}

# Find the optimal permutation of estimated clusters.
assignment = solve_LSAP(cost)

# Minimum total squared error after label matching.
min_sse = sum(
  cost[cbind(1:KC, assignment)]
)

# RMSE of the estimated cluster means.
MSEgwhe = c(
  MSEgwhe,
  sqrt(min_sse / (KC * p))
)


### ============================================================
### Covariance estimation error after iterative reassignment
### ============================================================

# Reconstruct the true covariance matrix.
sig = matrix(0.3, p, p) + 0.7 * diag(p)

# Compute RMSE using only the unique elements of the symmetric
# covariance matrix.
MSEhesig = c(
  MSEhesig,
  sqrt(
    sum(
      (
        Csig[upper.tri(Csig, diag = TRUE)] -
          jj^2 / 100 *
          sig[upper.tri(sig, diag = TRUE)]
      )^2
    ) /
      (p * (p + 1) / 2)
  )
)


### ============================================================
### Store ARI
### ============================================================

rii2 = c(rii2, ri2)