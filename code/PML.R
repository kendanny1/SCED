
rm(list=ls())

library(clue)

# ============================================================
# Basic settings
# ============================================================

nl=250
nu=0L                    # Number of unlabeled observations
nd=nl+nu                 # Total number of observations
KC=7                     # Number of clusters
p=15L                    # Number of attributes (dimensions)


# ============================================================
# Quartic kernel function
# ============================================================
# This is the quartic (biweight) kernel:
#
#   K(u) = (15/16)(1-u^2)^2,  |u| <= 1
#          0,                  otherwise
#
# It is used for kernel density estimation of the density
# generator.

K=function(u){
  out=15/16*(1-u^2)^2*(abs(u)<=1)
  return(out)
}


# ============================================================
# Helper functions for observations and matrices
# ============================================================

# Return column sums.
# If x is a vector of length p, return x itself.
# Otherwise, calculate the column sums of x.

cols2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=colSums(x)
  }
  return(out)
}


# Return the sum of all elements.
# If x is a vector of length p, return the sum of x.
# Otherwise, calculate the column sums of x.

cols3=function(x){
  if(length(x)==p){
    out=sum(x)
  }else{
    out=colSums(x)
  }
  return(out)
}


# Construct an outer-product matrix.
# If x is a vector of length p, return a p x p zero matrix.
# Otherwise, calculate x %*% t(x).

mo=function(x){
  if(length(x)==p){
    out=matrix(0,p,p)
  }else{
    out=x%*%t(x)
  }
  return(out)
}


# Transpose a matrix.
# If x is a vector of length p, return x itself.

t2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=t(x)
  }
  return(out)
}


# ============================================================
# Adjusted Rand Index (ARI)
# ============================================================
# This function calculates the Adjusted Rand Index between
# two cluster assignments.

adjustedRandIndex=function (x, y) 
{
  x <- as.vector(x)
  y <- as.vector(y)
  
  if (length(x) != length(y)) 
    stop("arguments must be vectors of the same length")
  
  tab <- table(x, y)
  
  if (all(dim(tab) == c(1, 1))) 
    return(1)
  
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
    ((a + b + a + c)/2 -
       (a + b) * (a + c)/(a + b + c + d))
  
  return(ARI)
}


# ============================================================
# Initialize simulation / experiment settings
# ============================================================

k=1
k2=6

# Initialize vectors used to store simulation results.
lfne=rie_1212=MSEkmsp_1212=MSEgw_1212=MSEsigel_1212=lfne=Clee2=c()

j1=1
KC=7


# ============================================================
# Load simulation results and data
# ============================================================

load(paste(
  "/NA3/kendanny1/sep/normalp15/ss/hsm_",
  KC,"_",k2,"_",nl,"_",500,".RData",
  sep=""
))

load(paste(
  "/NA3/kendanny1/sep/normalp15/data_",
  k2,".RData",
  sep=""
))

KC=7


# ============================================================
# Extract labeled observations and class labels
# ============================================================

# Extract the class labels for the j1-th simulation.
Cl=cxse[j1,]

# Extract the corresponding labeled observations.
Xl=Xlg[(nl*(j1-1)+1):(nl*j1),]


# ============================================================
# Sort observations according to their class labels
# ============================================================

D=cbind(Xl,Cl)

# Sort the data by the class label in ascending order.
D=D[order(D[,p+1],decreasing=FALSE),]


# Construct a class-membership indicator matrix.
pl=outer(D[,p+1],1:KC,FUN="==")
dim(pl)=c(nl,KC)

# Number of observations in each class.
nc=colSums(pl)

# Cumulative class sizes.
# Nc is used to locate the observations belonging to
# each class in the sorted data matrix.
Nc=rep(0,KC+1)
Nc[2:(KC+1)]=cumsum(nc)


# ============================================================
# Initial estimates of the class-specific means
# ============================================================

Cmu=c()

for (i in 1:KC){
  
  # Extract observations belonging to class i.
  Dc=D[D[,p+1]==i,1:p]
  
  # Calculate the sample mean vector for class i.
  Cmu=c(Cmu,cols2(Dc)/nc[i])
}


# ============================================================
# Initial estimate of the common covariance matrix
# ============================================================

Csig=matrix(0,p,p)

for(i in 1:KC){
  
  # Extract observations belonging to class i.
  Dc=D[D[,p+1]==i,1:p]
  
  # Calculate the class-specific sample mean.
  mc=cols2(Dc)/nc[i]
  
  # Accumulate the within-class scatter matrix.
  Csig=Csig+mo(t2(Dc)-mc)
}

# Estimate the common covariance matrix.
Csig=Csig/nl


# ============================================================
# Calculate initial Mahalanobis-type radial distances
# ============================================================
# The transformed quantity Y is used to estimate the
# density generator through kernel density estimation.

Y=c()

for(i in 1:KC){ 
  
  # Extract observations from class i.
  Dc=D[D[,p+1]==i,1:p]
  
  # Calculate the squared Mahalanobis-type distance.
  ylc=cols3(
    (solve(Csig)%*%
       (t2(Dc)-Cmu[((i-1)*p+1):(i*p)])) *
      (t2(Dc)-Cmu[((i-1)*p+1):(i*p)])
  )
  
  # Apply the transformation used for the elliptical
  # distributional density generator.
  Y=c(
    Y,
    (ylc^(p/2)+1)^(2/p)-1
  )
}


# ============================================================
# AMISE-based bandwidth selection
# ============================================================

# Create a grid for evaluating the estimated density.
y_seq=seq(min(Y),max(Y),length=5000)


# ------------------------------------------------------------
# AMISE objective function
# ------------------------------------------------------------

l=function(h){
  
  # Pairwise kernel evaluations between observations.
  K_1=
    K(outer(Y,Y,FUN="-")/h)/h+
    K(outer(Y,Y,FUN="+")/h)/h
  
  dim(K_1)=c(nl,nl)
  
  # Remove diagonal terms for the leave-one-out estimate.
  diag(K_1)=0
  
  # Integrated term used in the AMISE criterion.
  A_2=2*sum(K_1)/(nl*(nl-1))
  
  
  # Kernel density estimate evaluated on y_seq.
  K_2=
    K(outer(y_seq,Y,FUN="-")/h)/h+
    K(outer(y_seq,Y,FUN="+")/h)/h
  
  f_hat_h_y_seq=rowSums(K_2)/nl
  
  
  # Leave-one-out correction.
  K_3=
    nl/(nl-1)*f_hat_h_y_seq-
    K_2/(nl-1)
  
  
  # Approximate the integrated squared density
  # and subtract the pairwise kernel term.
  out=
    sum(K_3^2)*
    (max(Y)-min(Y))/(nl*5000)-
    A_2
  
  return(out)
}


# Obtain the initial bandwidth using numerical optimization.
h_hat0=nlminb(5,l)$par

# Apply the dimension-dependent bandwidth adjustment.
h_hat0=h_hat0*nl^(3/80)


# ============================================================
# Initial estimate of the density generator
# ============================================================

g_hat0=function(y){
  
  # Evaluate the kernel density estimate at y.
  K_4=
    K(
      outer(
        (y^(p/2)+1)^(2/p)-1,
        Y,
        FUN="-"
      )/h_hat0
    )/h_hat0+
    K(
      outer(
        (y^(p/2)+1)^(2/p)-1,
        Y,
        FUN="+"
      )/h_hat0
    )/h_hat0
  
  # Transform the kernel estimate into the estimated
  # density generator.
  out=
    rowSums(K_4)*
    (y^(p/2)+1)^(2/p-1)*
    gamma(p/2)/
    (nl*pi^(p/2))
  
  out
}


# ============================================================
# Initialize parameters for pseudo-likelihood estimation
# ============================================================

# Cholesky decomposition of the initial covariance matrix.
cho=chol(Csig)


# Parameter vector:
#   1. KC mean vectors
#   2. Upper-triangular Cholesky parameters
#   3. Class mixing proportions
ini=c(
  Cmu,
  cho[upper.tri(cho,diag=TRUE)],
  nc/nl
)


# ============================================================
# Groupwise pseudo-likelihood
# ============================================================
# The following objective function estimates the common
# scatter / covariance structure while keeping the class
# means fixed.

Ls=function(th){
  
  # Re-estimate the density generator using the current
  # transformed observations.
  g_hat0=function(y){
    
    K_4=
      K(
        outer(
          (y^(p/2)+1)^(2/p)-1,
          Y,
          FUN="-"
        )/h_hat0
      )/h_hat0
    
    out=
      rowSums(K_4)*
      (y^(p/2)+1)^(2/p-1)*
      gamma(p/2)/
      (nl*pi^(p/2))
    
    out
  }
  
  
  # Copy the current parameter vector.
  inir=ini
  
  # Replace the Cholesky parameters by the parameters
  # currently being optimized.
  inir[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]=th
  
  
  # Update the transformed radial distances Y for all classes.
  for(i in 1:KC){
    
    # Extract the mean vector of class i.
    mu=inir[((i-1)*p+1):(i*p)]
    
    # Extract the Cholesky parameters.
    siv=inir[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    # Reconstruct the upper-triangular Cholesky matrix.
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    # Construct the inverse transformation matrix.
    U=t(solve(V))
    
    # Extract observations belonging to class i.
    Dc=D[D[,p+1]==i,1:p]
    
    # Calculate transformed radial distances.
    Y[
      (Nc[i]+1):(Nc[i+1])
    ]=
      (
        cols3(
          (U%*%(t2(Dc)-mu))^2
        )^(p/2)+1
      )^(2/p)-1
  }
  
  
  # Initialize the total pseudo-log-likelihood.
  lfyu=0
  
  
  # ----------------------------------------------------------
  # Contributions from classes 1 through KC-1
  # ----------------------------------------------------------
  
  for(i in 1:(KC-1)){
    
    # Class-specific mean.
    mu=ini[((i-1)*p+1):(i*p)]
    
    # Current scatter parameters.
    siv=th
    
    # Reconstruct the Cholesky matrix.
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    # Inverse transformation matrix.
    U=t(solve(V))
    
    # Extract class-i observations.
    Dc=D[D[,p+1]==i,1:p]
    
    # Squared radial distances before the final transformation.
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    # Pseudo-log-likelihood contribution.
    lf=
      sum(log(g_hat0(Yl)))+
      nc[i]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
  }
  
  
  # ----------------------------------------------------------
  # Contribution from the final class
  # ----------------------------------------------------------
  
  mu=ini[((KC-1)*p+1):(KC*p)]
  siv=th
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==KC,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[KC]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  
  # Return the negative average pseudo-log-likelihood
  # because nlminb performs minimization.
  -lfyu/nl
}


# ============================================================
# Class-specific mean optimization functions
# ============================================================
# L1-L7 optimize the mean vector of one class at a time.
# The common scatter matrix and all other class means are
# kept fixed.

L1=function(th){
  
  # Estimate the density generator at the current radial
  # distances.
  g_hat0=function(y){
    
    K_4=
      K(
        outer(
          (y^(p/2)+1)^(2/p)-1,
          Y,
          FUN="-"
        )/h_hat0
      )/h_hat0
    
    out=
      rowSums(K_4)*
      (y^(p/2)+1)^(2/p-1)*
      gamma(p/2)/
      (nl*pi^(p/2))
    
    out
  }
  
  
  inir=ini
  
  # Replace the mean vector of class 1.
  inir[c(1:p)]=th
  
  
  # Recalculate Y using the updated class mean.
  for(i in 1:KC){
    
    mu=inir[((i-1)*p+1):(i*p)]
    
    siv=inir[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==i,1:p]
    
    Y[
      (Nc[i]+1):(Nc[i+1])
    ]=
      (
        cols3(
          (U%*%(t2(Dc)-mu))^2
        )^(p/2)+1
      )^(2/p)-1
  }
  
  
  lfyu=0
  
  # Updated mean for class 1.
  mu=th[1:p]
  
  # Keep the common scatter matrix fixed.
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==1,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[1]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  
  # Contribution from the reference/final class KC.
  mu=ini[((KC-1)*p+1):(KC*p)]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==KC,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[KC]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  -lfyu/nl
}


# ============================================================
# Optimize the mean vector of class 2
# ============================================================

L2=function(th){
  
  g_hat0=function(y){
    K_4=
      K(
        outer(
          (y^(p/2)+1)^(2/p)-1,
          Y,
          FUN="-"
        )/h_hat0
      )/h_hat0
    
    out=
      rowSums(K_4)*
      (y^(p/2)+1)^(2/p-1)*
      gamma(p/2)/
      (nl*pi^(p/2))
    
    out
  }
  
  inir=ini
  
  # Replace the mean vector of class 2.
  inir[c((p+1):(2*p))]=th
  
  for(i in 1:KC){
    
    mu=inir[((i-1)*p+1):(i*p)]
    
    siv=inir[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==i,1:p]
    
    Y[
      (Nc[i]+1):(Nc[i+1])
    ]=
      (
        cols3(
          (U%*%(t2(Dc)-mu))^2
        )^(p/2)+1
      )^(2/p)-1
  }
  
  lfyu=0
  
  # Updated mean for class 2.
  mu=th[1:p]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==2,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[2]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  
  # Contribution from the reference/final class KC.
  mu=ini[((KC-1)*p+1):(KC*p)]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==KC,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[KC]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  -lfyu/nl
}


# ============================================================
# Optimize the mean vector of class 3
# ============================================================

L3=function(th){
  
  g_hat0=function(y){
    K_4=
      K(
        outer(
          (y^(p/2)+1)^(2/p)-1,
          Y,
          FUN="-"
        )/h_hat0
      )/h_hat0
    
    out=
      rowSums(K_4)*
      (y^(p/2)+1)^(2/p-1)*
      gamma(p/2)/
      (nl*pi^(p/2))
    
    out
  }
  
  inir=ini
  
  # Replace the mean vector of class 3.
  inir[c((2*p+1):(3*p))]=th
  
  for(i in 1:KC){
    
    mu=inir[((i-1)*p+1):(i*p)]
    
    siv=inir[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==i,1:p]
    
    Y[
      (Nc[i]+1):(Nc[i+1])
    ]=
      (
        cols3(
          (U%*%(t2(Dc)-mu))^2
        )^(p/2)+1
      )^(2/p)-1
  }
  
  lfyu=0
  
  # Updated mean for class 3.
  mu=th[1:p]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==3,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[3]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  
  # Contribution from the reference/final class KC.
  mu=ini[((KC-1)*p+1):(KC*p)]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==KC,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[KC]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  -lfyu/nl
}


# ============================================================
# Optimize the mean vector of class 4
# ============================================================

L4=function(th){
  
  g_hat0=function(y){
    K_4=
      K(
        outer(
          (y^(p/2)+1)^(2/p)-1,
          Y,
          FUN="-"
        )/h_hat0
      )/h_hat0
    
    out=
      rowSums(K_4)*
      (y^(p/2)+1)^(2/p-1)*
      gamma(p/2)/
      (nl*pi^(p/2))
    
    out
  }
  
  inir=ini
  
  # Replace the mean vector of class 4.
  inir[c((3*p+1):(4*p))]=th
  
  for(i in 1:KC){
    
    mu=inir[((i-1)*p+1):(i*p)]
    
    siv=inir[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==i,1:p]
    
    Y[
      (Nc[i]+1):(Nc[i+1])
    ]=
      (
        cols3(
          (U%*%(t2(Dc)-mu))^2
        )^(p/2)+1
      )^(2/p)-1
  }
  
  lfyu=0
  
  # Updated mean for class 4.
  mu=th[1:p]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==4,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[4]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  
  # Contribution from the reference/final class KC.
  mu=ini[((KC-1)*p+1):(KC*p)]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==KC,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[KC]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  -lfyu/nl
}


# ============================================================
# Optimize the mean vector of class 5
# ============================================================

L5=function(th){
  
  g_hat0=function(y){
    K_4=
      K(
        outer(
          (y^(p/2)+1)^(2/p)-1,
          Y,
          FUN="-"
        )/h_hat0
      )/h_hat0
    
    out=
      rowSums(K_4)*
      (y^(p/2)+1)^(2/p-1)*
      gamma(p/2)/
      (nl*pi^(p/2))
    
    out
  }
  
  inir=ini
  
  # Replace the mean vector of class 5.
  inir[c((4*p+1):(5*p))]=th
  
  for(i in 1:KC){
    
    mu=inir[((i-1)*p+1):(i*p)]
    
    siv=inir[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==i,1:p]
    
    Y[
      (Nc[i]+1):(Nc[i+1])
    ]=
      (
        cols3(
          (U%*%(t2(Dc)-mu))^2
        )^(p/2)+1
      )^(2/p)-1
  }
  
  lfyu=0
  
  # Updated mean for class 5.
  mu=th[1:p]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==5,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[5]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  
  # Contribution from the reference/final class KC.
  mu=ini[((KC-1)*p+1):(KC*p)]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==KC,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[KC]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  -lfyu/nl
}


# ============================================================
# Optimize the mean vector of class 6
# ============================================================

L6=function(th){
  
  g_hat0=function(y){
    K_4=
      K(
        outer(
          (y^(p/2)+1)^(2/p)-1,
          Y,
          FUN="-"
        )/h_hat0
      )/h_hat0
    
    out=
      rowSums(K_4)*
      (y^(p/2)+1)^(2/p-1)*
      gamma(p/2)/
      (nl*pi^(p/2))
    
    out
  }
  
  inir=ini
  
  # Replace the mean vector of class 6.
  inir[c((5*p+1):(6*p))]=th
  
  for(i in 1:KC){
    
    mu=inir[((i-1)*p+1):(i*p)]
    
    siv=inir[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==i,1:p]
    
    Y[
      (Nc[i]+1):(Nc[i+1])
    ]=
      (
        cols3(
          (U%*%(t2(Dc)-mu))^2
        )^(p/2)+1
      )^(2/p)-1
  }
  
  lfyu=0
  
  # Updated mean for class 6.
  mu=th[1:p]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==6,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[6]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  
  # Contribution from the reference/final class KC.
  mu=ini[((KC-1)*p+1):(KC*p)]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==KC,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[KC]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  -lfyu/nl
}


# ============================================================
# Optimize the mean vector of class 7
# ============================================================

L7=function(th){
  
  g_hat0=function(y){
    K_4=
      K(
        outer(
          (y^(p/2)+1)^(2/p)-1,
          Y,
          FUN="-"
        )/h_hat0
      )/h_hat0
    
    out=
      rowSums(K_4)*
      (y^(p/2)+1)^(2/p-1)*
      gamma(p/2)/
      (nl*pi^(p/2))
    
    out
  }
  
  inir=ini
  
  # Replace the mean vector of class 7.
  inir[(6*p+1):(7*p)]=th
  
  for(i in 1:KC){
    
    mu=inir[((i-1)*p+1):(i*p)]
    
    siv=inir[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==i,1:p]
    
    Y[
      (Nc[i]+1):(Nc[i+1])
    ]=
      (
        cols3(
          (U%*%(t2(Dc)-mu))^2
        )^(p/2)+1
      )^(2/p)-1
  }
  
  lfyu=0
  
  # Updated mean for class 7.
  mu=th[1:p]
  
  siv=ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]
  
  V=matrix(0,p,p)
  V[upper.tri(V,diag=TRUE)]=siv
  
  U=t(solve(V))
  
  Dc=D[D[,p+1]==KC,1:p]
  
  Yl=cols3(
    (U%*%(t2(Dc)-mu))^2
  )
  
  lf=
    sum(log(g_hat0(Yl)))+
    nc[KC]*(sum(0.5*log(diag(U)^2)))
  
  lfyu=lfyu+lf
  
  -lfyu/nl
}


# ============================================================
# Iterative estimation procedure
# ============================================================
# The parameters are estimated by alternating optimization:
#
#   Step 1: Update the common scatter matrix.
#   Step 2: Update the mean vector of class 1.
#   Step 3: Update the mean vector of class 2.
#   ...
#   Step 7: Update the mean vector of class 7.
#
# The class proportions are then updated so that they sum to 1.
#
# The procedure continues until the change in the objective
# function is sufficiently small.

lfo=0
lfn=10

while(abs(lfo-lfn)>0.01){
  
  # Store the previous objective value.
  lfo=lfn
  
  
  # ----------------------------------------------------------
  # Step 1: Estimate the common scatter matrix
  # ----------------------------------------------------------
  
  ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]=
    nlminb(
      ini[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ],
      Ls,
      control=list(
        trace=1,
        rel.tol=1e-3,
        x.tol=5e-4,
        iter.max=7
      )
    )$par
  
  
  # Update the final class proportion so that all
  # class proportions sum to one.
  ini[
    p*KC+p*(p+1)/2+KC
  ]=
    1-
    sum(
      ini[
        (p*KC+p*(p+1)/2+1):
          (p*KC+p*(p+1)/2+KC-1)
      ]
    )
  
  
  # ----------------------------------------------------------
  # Step 2: Update class-1 mean
  # ----------------------------------------------------------
  
  ini[c(1:p)]=
    nlminb(
      ini[c(1:p)],
      L1,
      control=list(
        trace=1,
        rel.tol=1e-3,
        x.tol=5e-4,
        iter.max=7
      )
    )$par
  
  ini[
    p*KC+p*(p+1)/2+KC
  ]=
    1-
    sum(
      ini[
        (p*KC+p*(p+1)/2+1):
          (p*KC+p*(p+1)/2+KC-1)
      ]
    )
  
  
  # ----------------------------------------------------------
  # Step 3: Update class-2 mean
  # ----------------------------------------------------------
  
  ini[c((p+1):(2*p))]=
    nlminb(
      ini[c((p+1):(2*p))],
      L2,
      control=list(
        trace=1,
        rel.tol=1e-3,
        x.tol=5e-4,
        iter.max=7
      )
    )$par
  
  ini[
    p*KC+p*(p+1)/2+KC
  ]=
    1-
    sum(
      ini[
        (p*KC+p*(p+1)/2+1):
          (p*KC+p*(p+1)/2+KC-1)
      ]
    )
  
  
  # ----------------------------------------------------------
  # Step 4: Update class-3 mean
  # ----------------------------------------------------------
  
  ini[c((2*p+1):(3*p))]=
    nlminb(
      ini[c((2*p+1):(3*p))],
      L3,
      control=list(
        trace=1,
        rel.tol=1e-3,
        x.tol=5e-4,
        iter.max=7
      )
    )$par
  
  ini[
    p*KC+p*(p+1)/2+KC
  ]=
    1-
    sum(
      ini[
        (p*KC+p*(p+1)/2+1):
          (p*KC+p*(p+1)/2+KC-1)
      ]
    )
  
  
  # ----------------------------------------------------------
  # Step 5: Update class-4 mean
  # ----------------------------------------------------------
  
  ini[c((3*p+1):(4*p))]=
    nlminb(
      ini[c((3*p+1):(4*p))],
      L4,
      control=list(
        trace=1,
        rel.tol=1e-3,
        x.tol=5e-4,
        iter.max=7
      )
    )$par
  
  ini[
    p*KC+p*(p+1)/2+KC
  ]=
    1-
    sum(
      ini[
        (p*KC+p*(p+1)/2+1):
          (p*KC+p*(p+1)/2+KC-1)
      ]
    )
  
  
  # ----------------------------------------------------------
  # Step 6: Update class-5 mean
  # ----------------------------------------------------------
  
  ini[(4*p+1):(5*p)]=
    nlminb(
      ini[(4*p+1):(5*p)],
      L5,
      control=list(
        trace=1,
        rel.tol=1e-3,
        x.tol=5e-4,
        iter.max=7
      )
    )$par
  
  ini[
    p*KC+p*(p+1)/2+KC
  ]=
    1-
    sum(
      ini[
        (p*KC+p*(p+1)/2+1):
          (p*KC+p*(p+1)/2+KC-1)
      ]
    )
  
  
  # ----------------------------------------------------------
  # Step 7: Update class-6 mean
  # ----------------------------------------------------------
  
  ini[(5*p+1):(6*p)]=
    nlminb(
      ini[(5*p+1):(6*p)],
      L6,
      control=list(
        trace=1,
        rel.tol=1e-3,
        x.tol=5e-4,
        iter.max=7
      )
    )$par
  
  ini[
    p*KC+p*(p+1)/2+KC
  ]=
    1-
    sum(
      ini[
        (p*KC+p*(p+1)/2+1):
          (p*KC+p*(p+1)/2+KC-1)
      ]
    )
  
  
  # ----------------------------------------------------------
  # Step 8: Update class-7 mean
  # ----------------------------------------------------------
  
  ini[(6*p+1):(7*p)]=
    nlminb(
      ini[(6*p+1):(7*p)],
      L7,
      control=list(
        trace=1,
        rel.tol=1e-3,
        x.tol=5e-4,
        iter.max=7
      )
    )$par
  
  ini[
    p*KC+p*(p+1)/2+KC
  ]=
    1-
    sum(
      ini[
        (p*KC+p*(p+1)/2+1):
          (p*KC+p*(p+1)/2+KC-1)
      ]
    )
  
  
  # Evaluate the updated objective function.
  lfn=
    Ls(
      ini[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ]
    )
}


# ============================================================
# Final bandwidth adjustment and objective evaluation
# ============================================================

# Reverse the temporary bandwidth adjustment.
h_hat0=h_hat0*nl^(-3/80)

# Evaluate the final pseudo-likelihood objective.
lfn=
  Ls(
    ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
  )

lfni=lfn


# ============================================================
# Pseudo-likelihood estimator
# ============================================================

# Reconstruct the estimated upper-triangular Cholesky matrix.
V2=matrix(0,p,p)

V2[
  upper.tri(V2,diag=TRUE)
]=
  ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]


# Estimated scatter matrix.
pssi=t(V2)%*%V2


# Determinant-related normalizing constant.
detp=det(pssi)^(-0.5)


# ============================================================
# Extract estimated class means
# ============================================================

pmu=Y2=c()

for(i in 1:KC){
  
  # Store the estimated mean vector for class i.
  pmu=c(
    pmu,
    ini[((i-1)*p+1):(i*p)]
  )
}


# ============================================================
# Calculate radial distances using the final estimates
# ============================================================

for(i in 1:KC){ 
  
  # Extract observations belonging to class i.
  Dc=D[D[,p+1]==i,1:p]
  
  # Calculate the final squared Mahalanobis-type distances.
  ylc=
    colSums(
      (
        solve(pssi)%*%
          (t2(Dc)-pmu[(p*(i-1)+1):(p*i)])
      )*
        (t2(Dc)-pmu[(p*(i-1)+1):(p*i)])
    )
  
  # Apply the radial transformation.
  Y2=c(
    Y2,
    (ylc^(p/2)+1)^(2/p)-1
  )
}


# ============================================================
# Classification
# ============================================================

ri3=0

# Yt stores the class-specific squared radial distances
# for every observation.
Yt=c()

for(i in 1:KC){
  
  # Calculate the squared Mahalanobis-type distance of every
  # observation in Xl from the estimated mean of class i.
  Yt=cbind(
    Yt,
    cols3(
      (
        solve(pssi)%*%
          (t(Xl)-pmu[(p*(i-1)+1):(p*i)])
      )*
        (t(Xl)-pmu[(p*(i-1)+1):(p*i)])
    )
  )
}


# ============================================================
# Calculate posterior-like class probabilities
# ============================================================

pp=matrix(0,nd,KC)

for (i in 1:KC){
  
  # Calculate the class-specific density contribution:
  #
  #   determinant term
  #   × estimated density generator
  #   × class mixing proportion
  
  pp[,i]=
    detp*
    g_hat0(Yt[,i])*
    ini[
      KC*p+p*(p+1)/2+i
    ]
}


# ============================================================
# Assign each observation to the most likely class
# ============================================================

Clee=c()

for(i in 1:nd){
  
  # Identify the class with the largest posterior-like
  # probability.
  a=which(
    pp[i,]==max(pp[i,])
  )
  
  
  # If several classes have exactly the same maximum
  # probability, randomly select one of them.
  if(length(a)>1){
    b=sample(
      a,
      1,
      replace=FALSE,
      prob=rep(1,length(a))
    )
  }
  
  
  # If there is a unique maximum, assign that class directly.
  if(length(a)==1){
    b=a
  }
  
  
  # Store the predicted class label.
  Clee=c(Clee,b)
}
l4c=0
while(sum(Clee==Cl)<nl){
  l4c=l4c+1
  #load data
  Cl=Clee
  Xl=Xlg[(nl*(j1-1)+1):(nl*j1),]
  # ============================================================
  # Sort observations according to their class labels
  # ============================================================
  
  D=cbind(Xl,Cl)
  
  # Sort the data by the class label in ascending order.
  D=D[order(D[,p+1],decreasing=FALSE),]
  
  
  # Construct a class-membership indicator matrix.
  pl=outer(D[,p+1],1:KC,FUN="==")
  dim(pl)=c(nl,KC)
  
  # Number of observations in each class.
  nc=colSums(pl)
  
  # Cumulative class sizes.
  # Nc is used to locate the observations belonging to
  # each class in the sorted data matrix.
  Nc=rep(0,KC+1)
  Nc[2:(KC+1)]=cumsum(nc)
  
  
  # ============================================================
  # Initial estimates of the class-specific means
  # ============================================================
  
  Cmu=c()
  
  for (i in 1:KC){
    
    # Extract observations belonging to class i.
    Dc=D[D[,p+1]==i,1:p]
    
    # Calculate the sample mean vector for class i.
    Cmu=c(Cmu,cols2(Dc)/nc[i])
  }
  
  
  # ============================================================
  # Initial estimate of the common covariance matrix
  # ============================================================
  
  Csig=matrix(0,p,p)
  
  for(i in 1:KC){
    
    # Extract observations belonging to class i.
    Dc=D[D[,p+1]==i,1:p]
    
    # Calculate the class-specific sample mean.
    mc=cols2(Dc)/nc[i]
    
    # Accumulate the within-class scatter matrix.
    Csig=Csig+mo(t2(Dc)-mc)
  }
  
  # Estimate the common covariance matrix.
  Csig=Csig/nl
  
  
  # ============================================================
  # Calculate initial Mahalanobis-type radial distances
  # ============================================================
  # The transformed quantity Y is used to estimate the
  # density generator through kernel density estimation.
  
  Y=c()
  
  for(i in 1:KC){ 
    
    # Extract observations from class i.
    Dc=D[D[,p+1]==i,1:p]
    
    # Calculate the squared Mahalanobis-type distance.
    ylc=cols3(
      (solve(Csig)%*%
         (t2(Dc)-Cmu[((i-1)*p+1):(i*p)])) *
        (t2(Dc)-Cmu[((i-1)*p+1):(i*p)])
    )
    
    # Apply the transformation used for the elliptical
    # distributional density generator.
    Y=c(
      Y,
      (ylc^(p/2)+1)^(2/p)-1
    )
  }
  
  
  # ============================================================
  # AMISE-based bandwidth selection
  # ============================================================
  
  # Create a grid for evaluating the estimated density.
  y_seq=seq(min(Y),max(Y),length=5000)
  
  
  # ------------------------------------------------------------
  # AMISE objective function
  # ------------------------------------------------------------
  
  l=function(h){
    
    # Pairwise kernel evaluations between observations.
    K_1=
      K(outer(Y,Y,FUN="-")/h)/h+
      K(outer(Y,Y,FUN="+")/h)/h
    
    dim(K_1)=c(nl,nl)
    
    # Remove diagonal terms for the leave-one-out estimate.
    diag(K_1)=0
    
    # Integrated term used in the AMISE criterion.
    A_2=2*sum(K_1)/(nl*(nl-1))
    
    
    # Kernel density estimate evaluated on y_seq.
    K_2=
      K(outer(y_seq,Y,FUN="-")/h)/h+
      K(outer(y_seq,Y,FUN="+")/h)/h
    
    f_hat_h_y_seq=rowSums(K_2)/nl
    
    
    # Leave-one-out correction.
    K_3=
      nl/(nl-1)*f_hat_h_y_seq-
      K_2/(nl-1)
    
    
    # Approximate the integrated squared density
    # and subtract the pairwise kernel term.
    out=
      sum(K_3^2)*
      (max(Y)-min(Y))/(nl*5000)-
      A_2
    
    return(out)
  }
  
  
  # Obtain the initial bandwidth using numerical optimization.
  h_hat0=nlminb(5,l)$par
  
  # Apply the dimension-dependent bandwidth adjustment.
  h_hat0=h_hat0*nl^(3/80)
  
  
  # ============================================================
  # Initial estimate of the density generator
  # ============================================================
  
  g_hat0=function(y){
    
    # Evaluate the kernel density estimate at y.
    K_4=
      K(
        outer(
          (y^(p/2)+1)^(2/p)-1,
          Y,
          FUN="-"
        )/h_hat0
      )/h_hat0+
      K(
        outer(
          (y^(p/2)+1)^(2/p)-1,
          Y,
          FUN="+"
        )/h_hat0
      )/h_hat0
    
    # Transform the kernel estimate into the estimated
    # density generator.
    out=
      rowSums(K_4)*
      (y^(p/2)+1)^(2/p-1)*
      gamma(p/2)/
      (nl*pi^(p/2))
    
    out
  }
  
  
  # ============================================================
  # Initialize parameters for pseudo-likelihood estimation
  # ============================================================
  
  # Cholesky decomposition of the initial covariance matrix.
  cho=chol(Csig)
  
  
  # Parameter vector:
  #   1. KC mean vectors
  #   2. Upper-triangular Cholesky parameters
  #   3. Class mixing proportions
  ini=c(
    Cmu,
    cho[upper.tri(cho,diag=TRUE)],
    nc/nl
  )
  
  
  # ============================================================
  # Groupwise pseudo-likelihood
  # ============================================================
  # The following objective function estimates the common
  # scatter / covariance structure while keeping the class
  # means fixed.
  
  Ls=function(th){
    
    # Re-estimate the density generator using the current
    # transformed observations.
    g_hat0=function(y){
      
      K_4=
        K(
          outer(
            (y^(p/2)+1)^(2/p)-1,
            Y,
            FUN="-"
          )/h_hat0
        )/h_hat0
      
      out=
        rowSums(K_4)*
        (y^(p/2)+1)^(2/p-1)*
        gamma(p/2)/
        (nl*pi^(p/2))
      
      out
    }
    
    
    # Copy the current parameter vector.
    inir=ini
    
    # Replace the Cholesky parameters by the parameters
    # currently being optimized.
    inir[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]=th
    
    
    # Update the transformed radial distances Y for all classes.
    for(i in 1:KC){
      
      # Extract the mean vector of class i.
      mu=inir[((i-1)*p+1):(i*p)]
      
      # Extract the Cholesky parameters.
      siv=inir[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ]
      
      # Reconstruct the upper-triangular Cholesky matrix.
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      
      # Construct the inverse transformation matrix.
      U=t(solve(V))
      
      # Extract observations belonging to class i.
      Dc=D[D[,p+1]==i,1:p]
      
      # Calculate transformed radial distances.
      Y[
        (Nc[i]+1):(Nc[i+1])
      ]=
        (
          cols3(
            (U%*%(t2(Dc)-mu))^2
          )^(p/2)+1
        )^(2/p)-1
    }
    
    
    # Initialize the total pseudo-log-likelihood.
    lfyu=0
    
    
    # ----------------------------------------------------------
    # Contributions from classes 1 through KC-1
    # ----------------------------------------------------------
    
    for(i in 1:(KC-1)){
      
      # Class-specific mean.
      mu=ini[((i-1)*p+1):(i*p)]
      
      # Current scatter parameters.
      siv=th
      
      # Reconstruct the Cholesky matrix.
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      
      # Inverse transformation matrix.
      U=t(solve(V))
      
      # Extract class-i observations.
      Dc=D[D[,p+1]==i,1:p]
      
      # Squared radial distances before the final transformation.
      Yl=cols3(
        (U%*%(t2(Dc)-mu))^2
      )
      
      # Pseudo-log-likelihood contribution.
      lf=
        sum(log(g_hat0(Yl)))+
        nc[i]*(sum(0.5*log(diag(U)^2)))
      
      lfyu=lfyu+lf
    }
    
    
    # ----------------------------------------------------------
    # Contribution from the final class
    # ----------------------------------------------------------
    
    mu=ini[((KC-1)*p+1):(KC*p)]
    siv=th
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==KC,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[KC]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    
    # Return the negative average pseudo-log-likelihood
    # because nlminb performs minimization.
    -lfyu/nl
  }
  
  
  # ============================================================
  # Class-specific mean optimization functions
  # ============================================================
  # L1-L7 optimize the mean vector of one class at a time.
  # The common scatter matrix and all other class means are
  # kept fixed.
  
  L1=function(th){
    
    # Estimate the density generator at the current radial
    # distances.
    g_hat0=function(y){
      
      K_4=
        K(
          outer(
            (y^(p/2)+1)^(2/p)-1,
            Y,
            FUN="-"
          )/h_hat0
        )/h_hat0
      
      out=
        rowSums(K_4)*
        (y^(p/2)+1)^(2/p-1)*
        gamma(p/2)/
        (nl*pi^(p/2))
      
      out
    }
    
    
    inir=ini
    
    # Replace the mean vector of class 1.
    inir[c(1:p)]=th
    
    
    # Recalculate Y using the updated class mean.
    for(i in 1:KC){
      
      mu=inir[((i-1)*p+1):(i*p)]
      
      siv=inir[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ]
      
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      
      U=t(solve(V))
      
      Dc=D[D[,p+1]==i,1:p]
      
      Y[
        (Nc[i]+1):(Nc[i+1])
      ]=
        (
          cols3(
            (U%*%(t2(Dc)-mu))^2
          )^(p/2)+1
        )^(2/p)-1
    }
    
    
    lfyu=0
    
    # Updated mean for class 1.
    mu=th[1:p]
    
    # Keep the common scatter matrix fixed.
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==1,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[1]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    
    # Contribution from the reference/final class KC.
    mu=ini[((KC-1)*p+1):(KC*p)]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==KC,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[KC]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    -lfyu/nl
  }
  
  
  # ============================================================
  # Optimize the mean vector of class 2
  # ============================================================
  
  L2=function(th){
    
    g_hat0=function(y){
      K_4=
        K(
          outer(
            (y^(p/2)+1)^(2/p)-1,
            Y,
            FUN="-"
          )/h_hat0
        )/h_hat0
      
      out=
        rowSums(K_4)*
        (y^(p/2)+1)^(2/p-1)*
        gamma(p/2)/
        (nl*pi^(p/2))
      
      out
    }
    
    inir=ini
    
    # Replace the mean vector of class 2.
    inir[c((p+1):(2*p))]=th
    
    for(i in 1:KC){
      
      mu=inir[((i-1)*p+1):(i*p)]
      
      siv=inir[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ]
      
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      
      U=t(solve(V))
      
      Dc=D[D[,p+1]==i,1:p]
      
      Y[
        (Nc[i]+1):(Nc[i+1])
      ]=
        (
          cols3(
            (U%*%(t2(Dc)-mu))^2
          )^(p/2)+1
        )^(2/p)-1
    }
    
    lfyu=0
    
    # Updated mean for class 2.
    mu=th[1:p]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==2,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[2]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    
    # Contribution from the reference/final class KC.
    mu=ini[((KC-1)*p+1):(KC*p)]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==KC,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[KC]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    -lfyu/nl
  }
  
  
  # ============================================================
  # Optimize the mean vector of class 3
  # ============================================================
  
  L3=function(th){
    
    g_hat0=function(y){
      K_4=
        K(
          outer(
            (y^(p/2)+1)^(2/p)-1,
            Y,
            FUN="-"
          )/h_hat0
        )/h_hat0
      
      out=
        rowSums(K_4)*
        (y^(p/2)+1)^(2/p-1)*
        gamma(p/2)/
        (nl*pi^(p/2))
      
      out
    }
    
    inir=ini
    
    # Replace the mean vector of class 3.
    inir[c((2*p+1):(3*p))]=th
    
    for(i in 1:KC){
      
      mu=inir[((i-1)*p+1):(i*p)]
      
      siv=inir[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ]
      
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      
      U=t(solve(V))
      
      Dc=D[D[,p+1]==i,1:p]
      
      Y[
        (Nc[i]+1):(Nc[i+1])
      ]=
        (
          cols3(
            (U%*%(t2(Dc)-mu))^2
          )^(p/2)+1
        )^(2/p)-1
    }
    
    lfyu=0
    
    # Updated mean for class 3.
    mu=th[1:p]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==3,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[3]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    
    # Contribution from the reference/final class KC.
    mu=ini[((KC-1)*p+1):(KC*p)]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==KC,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[KC]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    -lfyu/nl
  }
  
  
  # ============================================================
  # Optimize the mean vector of class 4
  # ============================================================
  
  L4=function(th){
    
    g_hat0=function(y){
      K_4=
        K(
          outer(
            (y^(p/2)+1)^(2/p)-1,
            Y,
            FUN="-"
          )/h_hat0
        )/h_hat0
      
      out=
        rowSums(K_4)*
        (y^(p/2)+1)^(2/p-1)*
        gamma(p/2)/
        (nl*pi^(p/2))
      
      out
    }
    
    inir=ini
    
    # Replace the mean vector of class 4.
    inir[c((3*p+1):(4*p))]=th
    
    for(i in 1:KC){
      
      mu=inir[((i-1)*p+1):(i*p)]
      
      siv=inir[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ]
      
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      
      U=t(solve(V))
      
      Dc=D[D[,p+1]==i,1:p]
      
      Y[
        (Nc[i]+1):(Nc[i+1])
      ]=
        (
          cols3(
            (U%*%(t2(Dc)-mu))^2
          )^(p/2)+1
        )^(2/p)-1
    }
    
    lfyu=0
    
    # Updated mean for class 4.
    mu=th[1:p]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==4,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[4]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    
    # Contribution from the reference/final class KC.
    mu=ini[((KC-1)*p+1):(KC*p)]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==KC,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[KC]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    -lfyu/nl
  }
  
  
  # ============================================================
  # Optimize the mean vector of class 5
  # ============================================================
  
  L5=function(th){
    
    g_hat0=function(y){
      K_4=
        K(
          outer(
            (y^(p/2)+1)^(2/p)-1,
            Y,
            FUN="-"
          )/h_hat0
        )/h_hat0
      
      out=
        rowSums(K_4)*
        (y^(p/2)+1)^(2/p-1)*
        gamma(p/2)/
        (nl*pi^(p/2))
      
      out
    }
    
    inir=ini
    
    # Replace the mean vector of class 5.
    inir[c((4*p+1):(5*p))]=th
    
    for(i in 1:KC){
      
      mu=inir[((i-1)*p+1):(i*p)]
      
      siv=inir[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ]
      
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      
      U=t(solve(V))
      
      Dc=D[D[,p+1]==i,1:p]
      
      Y[
        (Nc[i]+1):(Nc[i+1])
      ]=
        (
          cols3(
            (U%*%(t2(Dc)-mu))^2
          )^(p/2)+1
        )^(2/p)-1
    }
    
    lfyu=0
    
    # Updated mean for class 5.
    mu=th[1:p]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==5,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[5]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    
    # Contribution from the reference/final class KC.
    mu=ini[((KC-1)*p+1):(KC*p)]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==KC,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[KC]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    -lfyu/nl
  }
  
  
  # ============================================================
  # Optimize the mean vector of class 6
  # ============================================================
  
  L6=function(th){
    
    g_hat0=function(y){
      K_4=
        K(
          outer(
            (y^(p/2)+1)^(2/p)-1,
            Y,
            FUN="-"
          )/h_hat0
        )/h_hat0
      
      out=
        rowSums(K_4)*
        (y^(p/2)+1)^(2/p-1)*
        gamma(p/2)/
        (nl*pi^(p/2))
      
      out
    }
    
    inir=ini
    
    # Replace the mean vector of class 6.
    inir[c((5*p+1):(6*p))]=th
    
    for(i in 1:KC){
      
      mu=inir[((i-1)*p+1):(i*p)]
      
      siv=inir[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ]
      
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      
      U=t(solve(V))
      
      Dc=D[D[,p+1]==i,1:p]
      
      Y[
        (Nc[i]+1):(Nc[i+1])
      ]=
        (
          cols3(
            (U%*%(t2(Dc)-mu))^2
          )^(p/2)+1
        )^(2/p)-1
    }
    
    lfyu=0
    
    # Updated mean for class 6.
    mu=th[1:p]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==6,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[6]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    
    # Contribution from the reference/final class KC.
    mu=ini[((KC-1)*p+1):(KC*p)]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==KC,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[KC]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    -lfyu/nl
  }
  
  
  # ============================================================
  # Optimize the mean vector of class 7
  # ============================================================
  
  L7=function(th){
    
    g_hat0=function(y){
      K_4=
        K(
          outer(
            (y^(p/2)+1)^(2/p)-1,
            Y,
            FUN="-"
          )/h_hat0
        )/h_hat0
      
      out=
        rowSums(K_4)*
        (y^(p/2)+1)^(2/p-1)*
        gamma(p/2)/
        (nl*pi^(p/2))
      
      out
    }
    
    inir=ini
    
    # Replace the mean vector of class 7.
    inir[(6*p+1):(7*p)]=th
    
    for(i in 1:KC){
      
      mu=inir[((i-1)*p+1):(i*p)]
      
      siv=inir[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ]
      
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      
      U=t(solve(V))
      
      Dc=D[D[,p+1]==i,1:p]
      
      Y[
        (Nc[i]+1):(Nc[i+1])
      ]=
        (
          cols3(
            (U%*%(t2(Dc)-mu))^2
          )^(p/2)+1
        )^(2/p)-1
    }
    
    lfyu=0
    
    # Updated mean for class 7.
    mu=th[1:p]
    
    siv=ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
    
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    
    U=t(solve(V))
    
    Dc=D[D[,p+1]==KC,1:p]
    
    Yl=cols3(
      (U%*%(t2(Dc)-mu))^2
    )
    
    lf=
      sum(log(g_hat0(Yl)))+
      nc[KC]*(sum(0.5*log(diag(U)^2)))
    
    lfyu=lfyu+lf
    
    -lfyu/nl
  }
  
  
  # ============================================================
  # Iterative estimation procedure
  # ============================================================
  # The parameters are estimated by alternating optimization:
  #
  #   Step 1: Update the common scatter matrix.
  #   Step 2: Update the mean vector of class 1.
  #   Step 3: Update the mean vector of class 2.
  #   ...
  #   Step 7: Update the mean vector of class 7.
  #
  # The class proportions are then updated so that they sum to 1.
  #
  # The procedure continues until the change in the objective
  # function is sufficiently small.
  
  lfo=0
  lfn=10
  
  while(abs(lfo-lfn)>0.01){
    
    # Store the previous objective value.
    lfo=lfn
    
    
    # ----------------------------------------------------------
    # Step 1: Estimate the common scatter matrix
    # ----------------------------------------------------------
    
    ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]=
      nlminb(
        ini[
          (p*KC+1):
            (p*KC+p*(p+1)/2)
        ],
        Ls,
        control=list(
          trace=1,
          rel.tol=1e-3,
          x.tol=5e-4,
          iter.max=7
        )
      )$par
    
    
    # Update the final class proportion so that all
    # class proportions sum to one.
    ini[
      p*KC+p*(p+1)/2+KC
    ]=
      1-
      sum(
        ini[
          (p*KC+p*(p+1)/2+1):
            (p*KC+p*(p+1)/2+KC-1)
        ]
      )
    
    
    # ----------------------------------------------------------
    # Step 2: Update class-1 mean
    # ----------------------------------------------------------
    
    ini[c(1:p)]=
      nlminb(
        ini[c(1:p)],
        L1,
        control=list(
          trace=1,
          rel.tol=1e-3,
          x.tol=5e-4,
          iter.max=7
        )
      )$par
    
    ini[
      p*KC+p*(p+1)/2+KC
    ]=
      1-
      sum(
        ini[
          (p*KC+p*(p+1)/2+1):
            (p*KC+p*(p+1)/2+KC-1)
        ]
      )
    
    
    # ----------------------------------------------------------
    # Step 3: Update class-2 mean
    # ----------------------------------------------------------
    
    ini[c((p+1):(2*p))]=
      nlminb(
        ini[c((p+1):(2*p))],
        L2,
        control=list(
          trace=1,
          rel.tol=1e-3,
          x.tol=5e-4,
          iter.max=7
        )
      )$par
    
    ini[
      p*KC+p*(p+1)/2+KC
    ]=
      1-
      sum(
        ini[
          (p*KC+p*(p+1)/2+1):
            (p*KC+p*(p+1)/2+KC-1)
        ]
      )
    
    
    # ----------------------------------------------------------
    # Step 4: Update class-3 mean
    # ----------------------------------------------------------
    
    ini[c((2*p+1):(3*p))]=
      nlminb(
        ini[c((2*p+1):(3*p))],
        L3,
        control=list(
          trace=1,
          rel.tol=1e-3,
          x.tol=5e-4,
          iter.max=7
        )
      )$par
    
    ini[
      p*KC+p*(p+1)/2+KC
    ]=
      1-
      sum(
        ini[
          (p*KC+p*(p+1)/2+1):
            (p*KC+p*(p+1)/2+KC-1)
        ]
      )
    
    
    # ----------------------------------------------------------
    # Step 5: Update class-4 mean
    # ----------------------------------------------------------
    
    ini[c((3*p+1):(4*p))]=
      nlminb(
        ini[c((3*p+1):(4*p))],
        L4,
        control=list(
          trace=1,
          rel.tol=1e-3,
          x.tol=5e-4,
          iter.max=7
        )
      )$par
    
    ini[
      p*KC+p*(p+1)/2+KC
    ]=
      1-
      sum(
        ini[
          (p*KC+p*(p+1)/2+1):
            (p*KC+p*(p+1)/2+KC-1)
        ]
      )
    
    
    # ----------------------------------------------------------
    # Step 6: Update class-5 mean
    # ----------------------------------------------------------
    
    ini[(4*p+1):(5*p)]=
      nlminb(
        ini[(4*p+1):(5*p)],
        L5,
        control=list(
          trace=1,
          rel.tol=1e-3,
          x.tol=5e-4,
          iter.max=7
        )
      )$par
    
    ini[
      p*KC+p*(p+1)/2+KC
    ]=
      1-
      sum(
        ini[
          (p*KC+p*(p+1)/2+1):
            (p*KC+p*(p+1)/2+KC-1)
        ]
      )
    
    
    # ----------------------------------------------------------
    # Step 7: Update class-6 mean
    # ----------------------------------------------------------
    
    ini[(5*p+1):(6*p)]=
      nlminb(
        ini[(5*p+1):(6*p)],
        L6,
        control=list(
          trace=1,
          rel.tol=1e-3,
          x.tol=5e-4,
          iter.max=7
        )
      )$par
    
    ini[
      p*KC+p*(p+1)/2+KC
    ]=
      1-
      sum(
        ini[
          (p*KC+p*(p+1)/2+1):
            (p*KC+p*(p+1)/2+KC-1)
        ]
      )
    
    
    # ----------------------------------------------------------
    # Step 8: Update class-7 mean
    # ----------------------------------------------------------
    
    ini[(6*p+1):(7*p)]=
      nlminb(
        ini[(6*p+1):(7*p)],
        L7,
        control=list(
          trace=1,
          rel.tol=1e-3,
          x.tol=5e-4,
          iter.max=7
        )
      )$par
    
    ini[
      p*KC+p*(p+1)/2+KC
    ]=
      1-
      sum(
        ini[
          (p*KC+p*(p+1)/2+1):
            (p*KC+p*(p+1)/2+KC-1)
        ]
      )
    
    
    # Evaluate the updated objective function.
    lfn=
      Ls(
        ini[
          (p*KC+1):
            (p*KC+p*(p+1)/2)
        ]
      )
  }
  
  
  # ============================================================
  # Final bandwidth adjustment and objective evaluation
  # ============================================================
  
  # Reverse the temporary bandwidth adjustment.
  h_hat0=h_hat0*nl^(-3/80)
  
  # Evaluate the final pseudo-likelihood objective.
  lfn=
    Ls(
      ini[
        (p*KC+1):
          (p*KC+p*(p+1)/2)
      ]
    )
  
  lfni=lfn
  
  
  # ============================================================
  # Pseudo-likelihood estimator
  # ============================================================
  
  # Reconstruct the estimated upper-triangular Cholesky matrix.
  V2=matrix(0,p,p)
  
  V2[
    upper.tri(V2,diag=TRUE)
  ]=
    ini[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
  
  
  # Estimated scatter matrix.
  pssi=t(V2)%*%V2
  
  
  # Determinant-related normalizing constant.
  detp=det(pssi)^(-0.5)
  
  
  # ============================================================
  # Extract estimated class means
  # ============================================================
  
  pmu=Y2=c()
  
  for(i in 1:KC){
    
    # Store the estimated mean vector for class i.
    pmu=c(
      pmu,
      ini[((i-1)*p+1):(i*p)]
    )
  }
  
  
  # ============================================================
  # Calculate radial distances using the final estimates
  # ============================================================
  
  for(i in 1:KC){ 
    
    # Extract observations belonging to class i.
    Dc=D[D[,p+1]==i,1:p]
    
    # Calculate the final squared Mahalanobis-type distances.
    ylc=
      colSums(
        (
          solve(pssi)%*%
            (t2(Dc)-pmu[(p*(i-1)+1):(p*i)])
        )*
          (t2(Dc)-pmu[(p*(i-1)+1):(p*i)])
      )
    
    # Apply the radial transformation.
    Y2=c(
      Y2,
      (ylc^(p/2)+1)^(2/p)-1
    )
  }
  
  
  # ============================================================
  # Classification
  # ============================================================
  
  ri3=0
  
  # Yt stores the class-specific squared radial distances
  # for every observation.
  Yt=c()
  
  for(i in 1:KC){
    
    # Calculate the squared Mahalanobis-type distance of every
    # observation in Xl from the estimated mean of class i.
    Yt=cbind(
      Yt,
      cols3(
        (
          solve(pssi)%*%
            (t(Xl)-pmu[(p*(i-1)+1):(p*i)])
        )*
          (t(Xl)-pmu[(p*(i-1)+1):(p*i)])
      )
    )
  }
  
  
  # ============================================================
  # Calculate posterior-like class probabilities
  # ============================================================
  
  pp=matrix(0,nd,KC)
  
  for (i in 1:KC){
    
    # Calculate the class-specific density contribution:
    #
    #   determinant term
    #   × estimated density generator
    #   × class mixing proportion
    
    pp[,i]=
      detp*
      g_hat0(Yt[,i])*
      ini[
        KC*p+p*(p+1)/2+i
      ]
  }
  
  
  # ============================================================
  # Assign each observation to the most likely class
  # ============================================================
  
  Clee=c()
  
  for(i in 1:nd){
    
    # Identify the class with the largest posterior-like
    # probability.
    a=which(
      pp[i,]==max(pp[i,])
    )
    
    
    # If several classes have exactly the same maximum
    # probability, randomly select one of them.
    if(length(a)>1){
      b=sample(
        a,
        1,
        replace=FALSE,
        prob=rep(1,length(a))
      )
    }
    
    
    # If there is a unique maximum, assign that class directly.
    if(length(a)==1){
      b=a
    }
    
    
    # Store the predicted class label.
    Clee=c(Clee,b)
  }
  
  if(l4c>=7) break
}
if(lfni<lfn) Cl=cxse[j1,]
ri3=adjustedRandIndex(Clee,Clg[(nl*(j1-1)+1):(nl*j1)])  
rie_1212=c(rie_1212,ri3)
mu_true <- cbind(mu1, mu2, mu3, mu4, mu5, mu6, mu7)
mu_est <- matrix(pmu,p,7)

# cost[i,j] = mu_est 的第 i 個向量配對 mu_true 的第 j 個向量時的 SSE
cost <- matrix(0, 7, 7)

for (i in 1:7) {
  for (j in 1:7) {
    cost[i, j] <- sum((mu_est[, i] - mu_true[, j])^2)
  }
}
assignment <- solve_LSAP(cost)
min_sse <- sum(cost[cbind(1:7, assignment)])
MSEgw_1212=c(MSEgw_1212,sqrt(min_sse/(KC*p)))
sig=matrix(0.3,p,p)+0.7*diag(p)
r_seq=seq(0,((max(Y)+1)^(p/2)-1)^(2/p),length=5000)
A[j]=sum(r_seq^(p+1)*g_hat0(r_seq^2)*2*pi^(p/2)/gamma(p/2)*((max(Y)+1)^(p/2)-1)^(2/p)/5000/p)
MSEsigel_1212=c(MSEsigel_1212,sqrt(sum((A[j]*pssi[upper.tri(pssi,diag=TRUE)]-k2^2/100*sig[upper.tri(sig,diag = TRUE)])^2)/(p*(p+1)/2)))
Clee2=rbind(Clee2,Clee)
lfne=c(lfne,lfn)
print(j1)
filename=paste(KC,"_",nl,"_",k2,sep="")
save.image(paste("/NA3/kendanny1/sep/normalp15/ss/ellip/he/he_",filename,"_",j1,".RData",sep=""),version = 2) 

