# ============================================================
# Model and simulation parameters
# ============================================================

KC=7                       # Number of clusters
nl=250                     # Number of labeled observations
k2=6                       # Simulation / data setting
nu=0L                      # Number of unlabeled observations
nd=nl+nu                   # Total number of observations
p=15L                      # Number of attributes (dimensions)


# ============================================================
# Kernel functions
# ============================================================

# Quartic (biweight) kernel:
#
#   K(u) = (15/16)(1-u^2)^2,  |u| <= 1
#          0,                  otherwise
#
# This kernel is used for kernel density estimation.

K=function(u){
  out=15/16*(1-u^2)^2*(abs(u)<=1)
  return(out)
}


# Fourth-order kernel function.
# This function is defined for higher-order kernel estimation
# and can be used to reduce the asymptotic bias of the density
# estimator.

K4=function(u){
  (1-3*u^2)*(1-u^2)^2*(abs(u)<=1)*105/64
}


# ============================================================
# Helper functions
# ============================================================

# Return the input vector if x is a vector of length p.
# Otherwise, calculate the column sums of x.

cols2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=colSums(x)
  }
  return(out)
}


# Return the sum of a vector of length p.
# Otherwise, calculate the column sums of x.

cols3=function(x){
  if(length(x)==p){
    out=sum(x)
  }else{
    out=colSums(x)
  }
  return(out)
}


# Calculate the outer-product matrix x %*% t(x).
# If x is a vector of length p, return a p x p zero matrix.

mo=function(x){
  if(length(x)==p){
    out=matrix(0,p,p)
  }else{
    out=x%*%t(x)
  }
  return(out)
}


# Transpose a matrix.
# If x is already a vector of length p, return x unchanged.

t2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=t(x)
  }
  return(out)
}


# ============================================================
# Initialize result containers
# ============================================================

A=rep(0,1)

# Store the pseudo-likelihood objective values.
lf_0314=c()

# Initialize containers for different performance measures.
rie_1212=
  MSEkmsp_1212=
  MSEgw_1212=
  MSEsigel_1212=
  lfne=
  Clee2=c()

j1=1


# ============================================================
# Load simulation results and data
# ============================================================

KC=7

# Load the previously obtained classification results /
# simulation objects.
load(
  paste(
    "/NA3/kendanny1/sep/normalp15/ss/ellip/he/he_",
    KC,"_",nl,"_",k2,"_",k*tim,".RData",
    sep=""
  )
)

# Load the simulated data corresponding to setting k2.
load(
  paste(
    "/NA3/kendanny1/sep/normalp15/data_",
    k2,".RData",
    sep=""
  )
)

KC=7


# ============================================================
# Extract labeled observations and class labels
# ============================================================

# Extract the class labels corresponding to the current
# simulation replicate.
Cl=Clee2[(j1-1)%%tim+1,]

# Extract the labeled observations for replicate j1.
Xl=Xlg[
  (nl*(j1-1)+1):(nl*j1),
]


# ============================================================
# Sort observations according to their class labels
# ============================================================

# Combine the observations and their class labels.
D=cbind(Xl,Cl)

# Sort observations by the class label in ascending order.
D=D[
  order(D[,p+1],decreasing=FALSE),
]

# Construct the class-membership indicator matrix.
pl=outer(
  D[,p+1],
  1:KC,
  FUN="=="
)

dim(pl)=c(nl,KC)

# Number of observations in each class.
nc=colSums(pl)

# Cumulative class sizes.
# Nc[i] and Nc[i+1] identify the range of observations
# belonging to class i in the sorted data matrix.
Nc=rep(0,KC+1)
Nc[2:(KC+1)]=cumsum(nc)


# ============================================================
# Initial labeled estimates of the class means
# ============================================================

Cmu=c()

for (i in 1:KC){
  
  # Extract observations belonging to class i.
  Dc=D[D[,p+1]==i,1:p]
  
  # Calculate the sample mean vector for class i.
  Cmu=c(
    Cmu,
    cols2(Dc)/nc[i]
  )
}


# ============================================================
# Initial estimate of the common covariance matrix
# ============================================================

Csig=matrix(0,p,p)

for(i in 1:KC){
  
  # Extract observations from class i.
  Dc=D[D[,p+1]==i,1:p]
  
  # Calculate the sample mean vector of class i.
  mc=cols2(Dc)/nc[i]
  
  # Accumulate the within-class scatter matrix.
  Csig=
    Csig+
    mo(t2(Dc)-mc)
}

# Estimate the common covariance matrix by dividing the
# total within-class scatter by the total sample size.
Csig=Csig/nl


# ============================================================
# Calculate initial radial distances
# ============================================================
# Y contains transformed squared Mahalanobis-type distances.
# These values are used to estimate the density generator.

Y=c()

for(i in 1:KC){ 
  
  # Extract observations from class i.
  Dc=D[D[,p+1]==i,1:p]
  
  # Calculate the squared Mahalanobis-type distance.
  ylc=
    cols3(
      (
        solve(Csig)%*%
          (t(Dc)-Cmu[((i-1)*p+1):(i*p)])
      )*
        (
          t(Dc)-Cmu[((i-1)*p+1):(i*p)])
    )
  )

# Apply the radial transformation.
Y=c(
  Y,
  (ylc^(p/2)+1)^(2/p)-1
)
}


# ============================================================
# AMISE-based bandwidth selection
# ============================================================

# Construct a fine grid for evaluating the density estimate.
y_seq=
  seq(
    min(Y),
    max(Y),
    length=5000
  )


# ------------------------------------------------------------
# AMISE objective function
# ------------------------------------------------------------

l=function(h){
  
  # Pairwise kernel evaluations among the observations.
  K_1=
    K(
      outer(Y,Y,FUN="-")/h
    )/h+
    K(
      outer(Y,Y,FUN="+")/h
    )/h
  
  dim(K_1)=c(nl,nl)
  
  # Remove diagonal terms for the leave-one-out calculation.
  diag(K_1)=0
  
  # Kernel-based integrated term in the AMISE criterion.
  A_2=
    2*sum(K_1)/(nl*(nl-1))
  
  
  # Evaluate the kernel density estimate on y_seq.
  K_2=
    K(
      outer(y_seq,Y,FUN="-")/h
    )/h+
    K(
      outer(y_seq,Y,FUN="+")/h
    )/h
  
  f_hat_h_y_seq=
    rowSums(K_2)/nl
  
  
  # Leave-one-out correction to the density estimate.
  K_3=
    nl/(nl-1)*f_hat_h_y_seq-
    K_2/(nl-1)
  
  
  # Approximate the AMISE criterion by numerical integration.
  out=
    sum(K_3^2)*
    (max(Y)-min(Y))/(nl*5000)-
    A_2
  
  return(out)
}


# Obtain the initial bandwidth by minimizing the AMISE
# objective function.
h_hat0=nlminb(5,l)$par

# Apply the dimension-dependent bandwidth adjustment.
h_hat0=h_hat0*nl^(3/80)


# ============================================================
# Initial estimate of the density generator
# ============================================================

g_hat0=function(y){
  
  # Evaluate the kernel density estimate at the transformed
  # radial distance y.
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
  
  # Convert the kernel estimate into an estimate of the
  # density generator for the elliptical distribution.
  out=
    rowSums(K_4)*
    (y^(p/2)+1)^(2/p-1)*
    gamma(p/2)/
    (nl*pi^(p/2))
  
  out
}


# ============================================================
# Initial parameter vector
# ============================================================

# Obtain the Cholesky factor of the initial covariance matrix.
cho=chol(Csig)

# The parameter vector contains:
#
#   1. KC class-specific mean vectors
#   2. Upper-triangular elements of the Cholesky factor
#   3. KC class mixing proportions
#
ini=
  c(
    Cmu,
    cho[upper.tri(cho,diag=TRUE)],
    nc/nl
  )


# ============================================================
# Full groupwise pseudo-likelihood objective function
# ============================================================
# The function L simultaneously estimates:
#
#   - all class-specific mean vectors,
#   - the common scatter matrix,
#   - the first KC-1 mixing proportions.
#
# The final mixing proportion is constrained to make all
# mixing proportions sum to one.

L=function(th){
  
  # ----------------------------------------------------------
  # Extract class-specific mean vectors
  # ----------------------------------------------------------
  
  mul=
    th[
      1:(p*KC)
    ]
  
  
  # ----------------------------------------------------------
  # Extract the upper-triangular Cholesky parameters
  # ----------------------------------------------------------
  
  sigl=
    th[
      (p*KC+1):
        (p*KC+p*(p+1)/2)
    ]
  
  
  # Reconstruct the upper-triangular Cholesky factor.
  sigt=matrix(0,p,p)
  
  sigt[
    upper.tri(sigt,diag=TRUE)
  ]=sigl
  
  
  # Construct the positive-definite scatter matrix.
  sigtt=
    t(sigt)%*%sigt
  
  
  # ----------------------------------------------------------
  # Extract the first KC-1 mixing proportions
  # ----------------------------------------------------------
  
  pic=
    th[
      (p*KC+p*(p+1)/2+1):
        (p*KC+p*(p+1)/2+KC-1)
    ]
  
  
  # The final mixing proportion is determined by the
  # constraint that all mixing proportions sum to one.
  pic=
    c(
      pic,
      1-sum(pic)
    )
  
  
  # Initialize the total log pseudo-likelihood.
  lf=0
  
  Y=c()
  
  
  # ==========================================================
  # Calculate class-specific radial distances
  # ==========================================================
  
  for(i in 1:KC){ 
    
    # Use all observations for classification likelihood.
    Dc=D[,1:p]
    
    # Calculate the squared Mahalanobis-type distance
    # relative to class i.
    Y=
      cbind(
        Y,
        cols2(
          (
            solve(sigtt)%*%
              (t(Dc)-mul[((i-1)*p+1):(i*p)])
          )*
            (
              t(Dc)-mul[((i-1)*p+1):(i*p)])
        )
      )
    )
  }
  
  
  # ==========================================================
  # Calculate the pseudo-log-likelihood
  # ==========================================================
  
  for(i in 1:nl){
    
    # Contribution of the determinant of the scatter matrix.
    #
    # The mixture density is obtained by summing the
    # class-specific density generators weighted by the
    # estimated mixing proportions.
    lf=
      lf-
      log(det(sigtt))/2+
      log(
        sum(
          g_hat0(Y[i,])*pic
        )
      )
  }
  
  
  # Return the negative average pseudo-log-likelihood,
  # because nlminb performs minimization.
  return(-lf/nl)
}


# ============================================================
# Estimate the model parameters
# ============================================================
# The last mixing proportion is excluded from the optimization
# because it is determined by the sum-to-one constraint.

ini[
  -(p*KC+p*(p+1)/2+KC)
]=
  nlminb(
    ini[
      -(p*KC+p*(p+1)/2+KC)
    ],
    L,
    control=list(
      trace=1,
      rel.tol=1e-5,
      x.tol=1e-5,
      iter.max=15
    )
  )$par


# Reverse the temporary bandwidth adjustment.
h_hat0=
  h_hat0*
  nl^(-3/80)


# Store the final pseudo-likelihood value.
lf_0314=
  c(
    lf_0314,
    L(
      ini[
        -(p*KC+p*(p+1)/2+KC)
      ]
    )
  )


# ============================================================
# Construct the final pseudo-likelihood estimator
# ============================================================

# Reconstruct the estimated upper-triangular Cholesky factor.
V2=matrix(0,p,p)

V2[
  upper.tri(V2,diag=TRUE)
]=
  ini[
    (p*KC+1):
      (p*KC+p*(p+1)/2)
  ]


# Construct the final estimated scatter matrix.
pssi=
  t(V2)%*%V2


# Determinant-related normalizing constant.
detp=
  det(pssi)^(-0.5)


# ============================================================
# Extract the estimated class means
# ============================================================

pmu=Y2=c()

for(i in 1:KC){
  
  # Append the estimated mean vector of class i.
  pmu=
    c(
      pmu,
      ini[
        ((i-1)*p+1):(i*p)
      ]
    )
}


# ============================================================
# Calculate final radial distances for labeled observations
# ============================================================

for(i in 1:KC){ 
  
  # Extract observations assigned to class i.
  Dc=
    D[
      D[,p+1]==i,
      1:p
    ]
  
  # Calculate squared Mahalanobis-type distances using
  # the final estimated scatter matrix and mean vector.
  ylc=
    cols3(
      (
        solve(pssi)%*%
          (
            t2(Dc)-
              pmu[(p*(i-1)+1):(p*i)]
          )
      )*
        (
          t2(Dc)-
            pmu[(p*(i-1)+1):(p*i)]
        )
    )
  
  # Apply the radial transformation.
  Y2=
    c(
      Y2,
      (ylc^(p/2)+1)^(2/p)-1
    )
}


# ============================================================
# Estimate class means from the final classification
# ============================================================

# Initialize the estimated class mean matrix.
# Each column corresponds to one cluster.
mu=
  matrix(
    0,
    p,
    KC
  )


# Initialize beta0.
# This appears to be used for a separate calculation based
# on the true class labels.
beta0=c()

for(i in 1:n){
  
  beta0=
    c(
      beta0,
      (Clg[(j1-1)*n+i]==1)*mu1+
        (Clg[(j1-1)*n+i]==2)*mu2
    )
}


# Calculate the empirical mean vector of each estimated cluster.
for(i in 1:KC) {
  
  # If cluster i contains exactly one observation,
  # directly assign that observation as the cluster mean.
  if(length(Xl[Clee==i,])==p){
    
    mu[,i]=
      Xl[Clee==i,]
    
  }else{
    
    # Otherwise, calculate the column-wise sample mean.
    mu[,i]=
      colSums(
        Xl[Clee==i,]
      )/
      dim(
        Xl[Clee==i,]
      )[1]
  }
}


# Store a copy of the estimated cluster means.
gamma=
  matrix(
    0,
    p,
    KC
  )

for(i in 1:KC)
  gamma[,i]=mu[,i]


# ============================================================
# True and estimated cluster means
# ============================================================

# Combine the true mean vectors of all seven clusters.
mu_true=
  cbind(
    mu1,
    mu2,
    mu3,
    mu4,
    mu5,
    mu6,
    mu7
  )


# Convert the estimated mean vector into a p x 7 matrix.
mu_est=
  matrix(
    pmu,
    p,
    7
  )


# ============================================================
# Cluster-label matching using the Hungarian algorithm
# ============================================================
# Cluster labels are arbitrary. Therefore, estimated clusters
# must be matched with the true clusters before calculating
# the mean squared error.
#
# The cost matrix contains the squared Euclidean distance
# between each estimated mean and each true mean.

cost=
  matrix(
    0,
    7,
    7
  )

for (i in 1:7) {
  for (j in 1:7) {
    
    cost[i, j]=
      sum(
        (
          mu_est[,i]-
            mu_true[,j]
        )^2
      )
  }
}


# Solve the minimum-cost assignment problem.
# The Hungarian algorithm finds the permutation of cluster
# labels that minimizes the total squared error.
assignment=
  solve_LSAP(cost)


# Calculate the minimum total sum of squared errors
# after optimal cluster-label matching.
min_sse=
  sum(
    cost[
      cbind(
        1:7,
        assignment
      )
    ]
  )


# Calculate the root mean squared error of the estimated
# cluster means, normalized by the total number of
# parameters KC * p.
MSEgw_1212=
  c(
    MSEgw_1212,
    sqrt(
      min_sse/(KC*p)
    )
  )


# ============================================================
# Evaluation of the scatter matrix estimation
# ============================================================

# Define the target scatter matrix used in the simulation.
sig=
  matrix(
    0.3,
    p,
    p
  )+
  0.7*diag(p)


# Construct a grid of radial distances for numerical
# integration.
r_seq=
  seq(
    0,
    ((max(Y)+1)^(p/2)-1)^(2/p),
    length=5000
  )


# Estimate the scale-related quantity A[j1] through
# numerical integration of the estimated density generator.
A[j1]=
  sum(
    r_seq^(p+1)*
      g_hat0(r_seq^2)*
      2*pi^(p/2)/
      gamma(p/2)*
      ((max(Y)+1)^(p/2)-1)^(2/p)/
      5000/
      p
  )


# Calculate the root mean squared error of the estimated
# upper-triangular scatter parameters.
#
# The estimated scatter parameters are compared with the
# corresponding entries of the target scatter matrix,
# adjusted by the estimated scale factor A[j1].
MSEsigel_1212=
  c(
    MSEsigel_1212,
    sqrt(
      sum(
        (
          A[j]*
            pssi[upper.tri(pssi,diag=TRUE)]-
            k2^2/100*
            sig[upper.tri(sig,diag=TRUE)]
        )^2
      )/
        (p*(p+1)/2)
    )
  )


# ============================================================
# Save the current simulation results
# ============================================================

# Construct the output file name from the number of clusters,
# simulation index, and data setting.
filename=
  paste(
    KC,
    "_",
    iii,
    "_",
    k2,
    sep=""
  )


# Save the complete R workspace for the current replicate.
save.image(
  paste(
    "/NA3/kendanny1/sep/normalp15/ss/ellip/em/",
    filename,
    "_",
    j1,
    ".RData",
    sep=""
  ),
  version=2
)