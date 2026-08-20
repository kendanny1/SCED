## ============================================================
## Clear workspace and load required package
## ============================================================

rm(list=ls())

library(expm)


## ============================================================
## Functions
## ============================================================

# Soft-thresholding function.
# This function is used in the ADMM update of eta.
#
# x : input value
# y : threshold
#
# The function shrinks |x| toward zero by y while preserving
# the sign of x. Values with |x| <= y are set to zero.
ST=function(x,y){
  (abs(x)-y)*(abs(x)-y>0)*sign(x)
}


# Calculate the column-wise sum of x.
#
# If x is already a p-dimensional vector, return x directly.
# Otherwise, calculate the column sums.
cols2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=colSums(x)
  }
  return(out)
}


## ============================================================
## Global parameters
## ============================================================

n=250                 # Number of observations
KC=7                  # Number of clusters used by the proposed method
KCn=7                 # Number of clusters in the initial k-means result
v=6                   # Simulation / dataset index
p=15                  # Number of variables
j=1                   # Current dataset index


## ============================================================
## Construct Cm matrix
## ============================================================

# Cm is an n x (n*p) indicator matrix.
#
# For observation i, the p elements corresponding to beta_i
# are all equal to 1.
#
# This matrix is used to convert the vectorized beta-gamma
# differences into observation-wise distances.

Cm=matrix(0,n,n*p)

for(i in 1:n)
  Cm[i,(p*(i-1)+1):(p*i)]=rep(1,p)


## ============================================================
## Initialize storage objects
## ============================================================

# Objects used to store the final results and performance
# measures for the current simulation setting.

c_1210=
  gamma_1210=
  MSEkmsp_1210=
  MSEhesp_1210=
  MSEkmsig_1210=
  MSEhesig_1210=
  MSEgw_1210=c()


## ============================================================
## Load clustering results and simulated data
## ============================================================

# Load the initial k-means clustering result.
# KCn determines which k-means result is loaded.

load(paste("km_",KCn,"_",v,"_",n,".RData",sep=""))

# Load the simulated dataset.
load(paste("data_",v,".RData",sep=""))


## ============================================================
## Initialize performance and result storage
## ============================================================

ril=c()       # Clustering agreement / evaluation measure
msel=c()      # Mean squared error of beta
S2=c()        # Final quadratic loss
SS=c()        # Objective function / loss values

cxx=c()       # Estimated cluster assignments
betat=c()     # Estimated observation-specific beta
gammat=c()    # Estimated cluster centers


## ============================================================
## Extract initial cluster assignment
## ============================================================

# Clek contains the clustering labels.
# Extract the labels corresponding to the j-th dataset.

cxt=Clek[((j-1)*n+1):(j*n)]


## ============================================================
## Estimate initial beta and gamma using k-means clusters
## ============================================================

beta=gamma=c()

# ADMM step-size related parameter
kappa=1/n

# Extract the observations corresponding to dataset j.
X=Xlg[((j-1)*n+1):(j*n),]


# Initial beta_i:
# Each observation is assigned the mean vector of its
# corresponding initial cluster.

for(i in 1:n)
  beta=c(
    beta,
    cols2(X[cxt==cxt[i],])/sum(cxt==cxt[i])
  )


# Initial cluster centers gamma_k:
# Calculate the mean vector of each initial cluster.

for(i in 1:KCn)
  gamma=cbind(
    gamma,
    cols2(X[cxt==i,])/sum(cxt==i)
  )


## ============================================================
## Estimate within-cluster covariance matrix
## ============================================================

Csig2=matrix(0,p,p)

for(i in 1:KCn){
  
  # Special case:
  # If cluster i contains only one observation,
  # X[cxt==i,] is a vector of length p.
  
  if(length(X[cxt==i,])==p){
    
    Csig2=
      Csig2+
      outer(
        X[cxt==i,]-gamma[,i],
        X[cxt==i,]-gamma[,i],
        FUN="*"
      )
    
  }else{
    
    # For clusters containing multiple observations,
    # calculate the within-cluster scatter matrix.
    
    Csig2=
      Csig2+
      (t(X[cxt==i,])-gamma[,i])%*%
      t(t(X[cxt==i,])-gamma[,i])
  }
}


# Normalize by the sample size.
Csig2=Csig2/n


## ============================================================
## Obtain initial beta when KCn < KC
## ============================================================

# If the initial k-means solution has fewer clusters than
# the desired number of clusters, load the k-means result
# corresponding to KC and use it to construct an initial beta.


if(KCn<KC){
  
  load(paste("km_",KC,"_",v,"_",n,".RData",sep=""))
  
  beta=gamma=c()
  
  X=Xlg[((j-1)*n+1):(j*n),]
  
  # Standardize the data using the estimated covariance matrix.
  X=X%*%solve(sqrtm(Csig2))
  
  # Initial beta based on the KC-cluster solution.
  for(i in 1:n)
    beta=c(
      beta,
      cols2(X[cxt==cxt[i],])/sum(cxt==cxt[i])
    )
  
  # Initial cluster centers.
  for(i in 1:KC)
    gamma=cbind(
      gamma,
      cols2(X[cxt==i,])/sum(cxt==i)
    )
  
  beta_ini=beta
}


## ============================================================
## Standardize the data
## ============================================================

X=Xlg[((j-1)*n+1):(j*n),]

# Whitening transformation using Csig2.
# After transformation, the estimated covariance matrix is
# approximately the identity matrix.

X=X%*%solve(sqrtm(Csig2))

# Covariance matrix in the standardized space.
Csig=diag(p)


## ============================================================
## Recalculate initial beta and gamma after standardization
## ============================================================

beta=gamma=c()


# Observation-specific initial beta.
for(i in 1:n)
  beta=c(
    beta,
    cols2(X[cxt==cxt[i],])/sum(cxt==cxt[i])
  )


# Initial cluster centers based on KCn clusters.
for(i in 1:KCn)
  gamma=cbind(
    gamma,
    cols2(X[cxt==i,])/sum(cxt==i)
  )


# Save initial cluster centers.
gamma_ini=gamma


## ============================================================
## Vectorized representation of beta, X and gamma
## ============================================================

beta2=beta

# Convert vectorized beta into a p x n matrix.
beta2=matrix(beta2,p,n)

# Vectorized representation of X.
X2=t(X)
dim(X2)=c(n*p,1)

# Repeat gamma for all n observations.
#
# gamma2 has dimension (n*p) x KCn.
# The block corresponding to observation i contains
# the same set of KCn cluster centers.

gamma2=do.call(
  rbind,
  replicate(n,gamma,simplify = FALSE)
)


## ============================================================
## Initial quadratic loss
## ============================================================

# Calculate the initial quadratic loss:
#
#   (1 / 2n) * || Sigma^{-1/2} (X - beta) ||_F^2
#
# Since Csig = I after standardization, this term represents
# the within-cluster fitting error.

Sn=
  sum(
    (sqrtm(solve(Csig))%*%(t(X)-beta2))^2
  )/(2*n)

SS=c(SS,Sn)


## ============================================================
## Store initial model performance
## ============================================================

# Initial clustering evaluation.
ril=c(ril,rii2[j])

# Initial mean squared error of beta.
msel=c(
  msel,
  sum((beta-beta0)^2)/(n*p)
)

# Store initial cluster assignment.
cxx=cbind(cxx,cxt)

# Store initial beta.
betat=cbind(betat,beta)

# Store initial gamma.
gammat=cbind(gammat,gamma)


## ============================================================
## Determine the initial beta for lambda selection
## ============================================================

# If the number of initial k-means clusters equals the target
# number of clusters, use the current beta as beta_ini.

if(KCn==KC)
  beta_ini=beta


## ============================================================
## Construct candidate lambda range
## ============================================================

# Calculate all absolute coordinate-wise differences between
# beta_ini and the initial cluster centers.
#
# lam_can contains the possible scale of lambda based on
# |beta_i - gamma_k|.

lam_can=abs(beta_ini-gamma2)


# Alternative construction of lam_can based on the total
# L1 distance between beta_i and gamma_k.
# This code is currently commented out.
#
# lam_can=c()
# for(i in 1:n)
#   for(il in 1:KCn)
#     lam_can=c(
#       lam_can,
#       sum(
#         abs(
#           beta_ini[((i-1)*p+1):(i*p)]-
#           gamma[,il]
#         )
#       )
#     )


# Lower bound of lambda.
#
# Exclude zero distances because a zero value would result
# in an unsuitable lower bound.

lam_1=
  n^0.01*
  min(lam_can[lam_can!=0])


# Upper bound of lambda.
lam_J=n^0.01*max(lam_can)


## ============================================================
## Main loop over candidate lambda values
## ============================================================

# Evaluate 10 candidate lambda values between lam_1 and lam_J.
#
# lambda is subsequently divided by n to account for the
# sample-size scaling used in the objective function.

for(lambda in seq(lam_1,lam_J,length=10)/n){
  
  
  ## ==========================================================
  ## Initialize beta, gamma and ADMM variables
  ## ==========================================================
  
  beta=gamma=c()
  
  kappa=1/n
  
  # Standardized data.
  X=
    Xlg[((j-1)*n+1):(j*n),]%
  %solve(sqrtm(Csig2))
  
  # Initialize dual variables for KCn clusters.
  nu=matrix(0,n*p,KCn)
  
  
  # Start the optimization from beta_ini.
  beta=beta_ini
  
  
  # Recalculate initial cluster centers.
  for(i in 1:KCn)
    gamma=cbind(
      gamma,
      cols2(X[cxt==i,])/sum(cxt==i)
    )
  
  
  ## ----------------------------------------------------------
  ## Prepare vectorized matrices
  ## ----------------------------------------------------------
  
  beta2=beta
  beta2=matrix(beta2,p,n)
  
  X2=t(X)
  dim(X2)=c(n*p,1)
  
  gamma2=
    do.call(
      rbind,
      replicate(n,gamma,simplify = FALSE)
    )
  
  
  ## ==========================================================
  ## Initialize eta and calculate initial objective
  ## ==========================================================
  
  # eta represents beta_i - gamma_k.
  eta=beta-gamma2
  
  
  # Calculate L1 distances between every observation-specific
  # beta_i and every cluster center gamma_k.
  
  A=Cm%*%abs(beta-gamma2)
  
  DD=c()
  
  
  # For each observation, retain the minimum distance to
  # any cluster center.
  
  for(i in 1:n)
    DD=c(DD,min(A[i,]))
  
  
  # Initial penalized objective:
  #
  #   quadratic loss + lambda * fusion penalty
  
  Sn=
    sum(
      (sqrtm(solve(Csig))%*%(t(X)-beta2))^2
    )/(2*n)+
    lambda*sum(DD)
  
  
  ## ==========================================================
  ## Initialize outer optimization loop
  ## ==========================================================
  
  S=100
  
  gamman=gamma
  betan=beta
  
  l1c=0
  
  # Store beta estimates from each iteration.
  beta3=c()
  
  # Store objective values from each iteration.
  SSS=c()
  
  # Store cluster assignments from each iteration.
  cxt1=c()
  
  # Store quadratic loss from each iteration.
  SS1=c()
  
  # Store gamma estimates from each iteration.
  gammab=c()
  
  
  ## ==========================================================
  ## Main optimization loop
  ## ==========================================================
  
  while(abs(Sn-S)>1e-10){
    
    l1c=l1c+1
    
    # Save the previous objective value.
    S=Sn
    
    
    ## --------------------------------------------------------
    ## Initialize dual variables
    ## --------------------------------------------------------
    
    nu=matrix(0,n*p,KCn)
    
    
    ## --------------------------------------------------------
    ## Current beta and gamma
    ## --------------------------------------------------------
    
    beta=betan[1:(n*p)]
    
    beta2=beta
    beta2=matrix(beta2,p,n)
    
    gamma=gamman
    
    
    # Difference between beta and gamma.
    eta=beta-gamma2
    
    
    ## ========================================================
    ## Calculate fusion penalty quantities
    ## ========================================================
    
    # A[i,k] is the L1 distance between beta_i and gamma_k.
    A=Cm%*%abs(beta-gamma2)
    
    
    # Sum distances to all cluster centers except k.
    A2=rowSums(A)-A
    
    
    # Sign of beta_i - gamma_k.
    B=
      ((beta-gamma2)>0)-
      ((beta-gamma2)<0)
    
    
    # Sum of signs over all clusters except k.
    B2=rowSums(B)-B
    
    
    ## --------------------------------------------------------
    ## Identify the closest alternative cluster
    ## --------------------------------------------------------
    
    C=matrix(0,n*p,KCn)
    
    for(i in 1:n)
      C[
        ((i-1)*p+1):(i*p),
        which.max(A2[i,])
      ]=1
    
    
    # Penalty-related correction term for the beta update.
    R=rowSums(C*B2)*lambda
    
    
    ## ========================================================
    ## Update beta
    ## ========================================================
    
    # Calculate the data-dependent component.
    X2s=solve(Csig)%*%t(X)
    X2s=X2s[1:(n*p)]
    
    
    # Construct the right-hand side of the beta update.
    bs=
      X2s/n+
      kappa*rowSums(gamma2)+
      kappa*rowSums(eta)-
      rowSums(nu)+
      R
    
    dim(bs)=c(p,n)
    
    
    # Closed-form update of beta.
    bt=
      solve(
        solve(Csig)/n+
          kappa*KCn*diag(p)
      )%*%bs
    
    betan=bt[1:(n*p)]
    
    
    ## ========================================================
    ## Update eta using soft-thresholding
    ## ========================================================
    
    etan=nun=matrix(0,p*n,KCn)
    
    
    # eta update for every cluster and every observation.
    #
    # ST() performs the soft-thresholding operation that
    # introduces the L1-type sparsity / fusion structure.
    
    for(i in 1:KCn)
      for(j2 in 1:n)
        etan[
          (1+(j2-1)*p):(j2*p),
          i
        ]=
      ST(
        betan[(1+(j2-1)*p):(j2*p)]-
          gamma2[(1+(j2-1)*p):(j2*p),i]+
          nu[(1+(j2-1)*p):(j2*p),i]/kappa,
        lambda/kappa
      )
    
    
    ## ========================================================
    ## Update dual variables
    ## ========================================================
    
    # ADMM dual-variable update.
    
    for(i in 1:KCn)
      nun[,i]=
      nu[,i]+
      kappa*
      (betan-gamma2[,i]-etan[,i])
    
    
    ## ========================================================
    ## Inner ADMM iteration
    ## ========================================================
    
    rr=c()       # Primal residual
    ss=c()       # Change in eta
    l2c=0        # Inner iteration counter
    
    
    # Continue until the ADMM constraint residual is sufficiently
    # small.
    
    while(
      sum(
        (betan[1:(n*p)]-gamma2-etan)^2
      )>n*1e-4
    ){
      
      l2c=l2c+1
      
      
      ## ------------------------------------------------------
      ## Update current beta, eta and dual variables
      ## ------------------------------------------------------
      
      beta=betan[1:(n*p)]
      
      beta2=beta
      beta2=matrix(beta2,p,n)
      
      eta=etan
      nu=nun
      
      
      ## ------------------------------------------------------
      ## Recalculate beta
      ## ------------------------------------------------------
      
      X2s=solve(Csig)%*%t(X)
      X2s=X2s[1:(n*p)]
      
      
      bs=
        X2s/n+
        kappa*rowSums(gamma2)+
        kappa*rowSums(eta)-
        rowSums(nu)+
        R
      
      dim(bs)=c(p,n)
      
      
      bt=
        solve(
          solve(Csig)/n+
            kappa*KCn*diag(p)
        )%*%bs
      
      betan=bt[1:(n*p)]
      
      
      ## ------------------------------------------------------
      ## Recalculate eta using soft-thresholding
      ## ------------------------------------------------------
      
      for(i in 1:KCn)
        for(j2 in 1:n)
          etan[
            (1+(j2-1)*p):(j2*p),
            i
          ]=
        ST(
          betan[(1+(j2-1)*p):(j2*p)]-
            gamma2[(1+(j2-1)*p):(j2*p),i]+
            nu[(1+(j2-1)*p):(j2*p),i]/kappa,
          lambda/kappa
        )
      
      
      ## ------------------------------------------------------
      ## Update dual variables
      ## ------------------------------------------------------
      
      for(i in 1:KCn)
        nun[,i]=
        nu[,i]+
        kappa*
        (betan-gamma2[,i]-etan[,i])
      
      
      ## ------------------------------------------------------
      ## Store convergence diagnostics
      ## ------------------------------------------------------
      
      rr=c(
        rr,
        sum(
          (betan[1:(n*p)]-gamma2-etan)^2
        )
      )
      
      ss=c(
        ss,
        sum((etan-eta)^2)
      )
      
      
      # Maximum of 30 inner ADMM iterations.
      if(l2c>=30) break
    }
    
    
    ## ========================================================
    ## Reassign observations to clusters
    ## ========================================================
    
    # Calculate the L1 distance between beta_i and each gamma_k.
    A=
      Cm%*%
      abs(betan[1:(n*p)]-gamma2)
    
    
    # Exclude the current cluster when searching for
    # the alternative cluster.
    A2=rowSums(A)-A
    
    
    # Sign information used in the penalty calculation.
    B=
      ((betan[1:(n*p)]-gamma2)>0)-
      ((betan[1:(n*p)]-gamma2)<0)
    
    B2=rowSums(B)-B
    
    
    # Determine the cluster associated with the smallest
    # effective distance.
    C2=c()
    
    for(i in 1:n)
      C2=c(
        C2,
        which.max(A2[i,])
      )
    
    
    ## ========================================================
    ## Update cluster centers gamma
    ## ========================================================
    
    betan2=betan
    dim(betan2)=c(p,n)
    
    gamman=matrix(0,p,KCn)
    
    
    for(i in 1:KCn){
      
      # If cluster i becomes empty, retain its initial
      # cluster center gamma_ini.
      if(sum(C2==i)==0){
        
        gamman[,i]=gamma_ini[,i]
        
        # If only one observation belongs to cluster i,
        # use that observation's beta directly.
      }else if(sum(C2==i)==1){
        
        gamman[,i]=betan2[,C2==i]
        
        # Otherwise, calculate the coordinate-wise median
        # of beta for all observations assigned to cluster i.
      }else{
        
        gamman[,i]=
          apply(
            betan2[,C2==i],
            MARGIN=1,
            FUN="median"
          )
      }
    }
    
    
    ## ========================================================
    ## Iterative refinement of gamma
    ## ========================================================
    
    # gammao stores the previous gamma.
    # The iteration continues until gamma converges.
    
    gammao=matrix(0,p,KCn)
    
    l3c=0
    
    
    while(
      sum(gammao-gamman)^2>
      KCn*1e-4
    ){
      
      l3c=l3c+1
      
      # Save current gamma before updating.
      gammao=gamman
      
      
      ## ------------------------------------------------------
      ## Repeat gamma for all observations
      ## ------------------------------------------------------
      
      gamma2=
        do.call(
          rbind,
          replicate(n,gammao,simplify = FALSE)
        )
      
      
      ## ------------------------------------------------------
      ## Recalculate observation-cluster distances
      ## ------------------------------------------------------
      
      A=
        Cm%*%
        abs(betan[1:(n*p)]-gamma2)
      
      A2=rowSums(A)-A
      
      
      ## ------------------------------------------------------
      ## Reassign observations
      ## ------------------------------------------------------
      
      C2=c()
      
      for(i in 1:n)
        C2=c(
          C2,
          which.max(A2[i,])
        )
      
      
      ## ------------------------------------------------------
      ## Update gamma using the new assignments
      ## ------------------------------------------------------
      
      betan2=betan
      dim(betan2)=c(p,n)
      
      gamman=matrix(0,p,KCn)
      
      
      for(i in 1:KCn){
        
        # Keep the previous gamma if the cluster is empty.
        if(sum(C2==i)==0){
          
          gamman[,i]=gammao[,i]
          
          # Use beta directly for a singleton cluster.
        }else if(sum(C2==i)==1){
          
          gamman[,i]=betan2[,C2==i]
          
          # Use coordinate-wise median for non-singleton clusters.
        }else{
          
          gamman[,i]=
            apply(
              betan2[,C2==i],
              MARGIN=1,
              FUN="median"
            )
        }
      }
      
      
      # Maximum of 30 gamma-refinement iterations.
      if(l3c>=30) break
    }
    
    
    ## ========================================================
    ## Recalculate objective function after gamma update
    ## ========================================================
    
    # Repeat the updated gamma for all observations.
    gamma2=
      do.call(
        rbind,
        replicate(n,gamman,simplify = FALSE)
      )
    
    
    # Calculate beta-gamma distances.
    A=
      Cm%*%
      abs(betan[1:(n*p)]-gamma2)
    
    DD=c()
    
    
    # Minimum distance from each observation to its closest
    # cluster center.
    for(i in 1:n)
      DD=c(
        DD,
        min(A[i,])
      )
    
    
    ## --------------------------------------------------------
    ## Calculate penalized objective function
    ## --------------------------------------------------------
    
    beta2=betan
    beta2=matrix(beta2,p,n)
    
    Sn=
      sum(
        (sqrtm(solve(Csig))%*%(t(X)-beta2))^2
      )/(2*n)+
      lambda*sum(DD)
    
    
    ## ========================================================
    ## Store intermediate results
    ## ========================================================
    
    # Store beta from the current iteration.
    beta3=cbind(beta3,betan)
    
    
    # Update gamma2 using the latest gamma.
    gamma2=
      do.call(
        rbind,
        replicate(n,gamman,simplify = FALSE)
      )
    
    
    # Recalculate distances.
    A=
      Cm%*%
      abs(betan[1:(n*p)]-gamma2)
    
    
    # Assign each observation to its closest cluster.
    cxt=
      apply(
        A,
        MARGIN=1,
        FUN="which.min"
      )
    
    
    # Store the total objective function.
    SSS=c(SSS,Sn)
    
    
    # Store the quadratic loss without the lambda penalty.
    SS1=c(
      SS1,
      sum(
        (sqrtm(solve(Csig))%*%(t(X)-beta2))^2
      )/(2*n)
    )
    
    
    # Store cluster assignments from the current iteration.
    cxt1=cbind(cxt1,cxt)
    
    
    # Store cluster centers.
    gammab=cbind(gammab,gamman)
    
    
    # Maximum of 30 outer optimization iterations.
    if(l1c>=30) break
  }
  
  
  ## ==========================================================
  ## Extract the final result for the current lambda
  ## ==========================================================
  
  # Final cluster assignment obtained at the last optimization
  # iteration.
  cxtn=cxt1[,l1c]
  
  
  # Corresponding quadratic loss.
  SS2=SS1[l1c]
  
  
  # Extract the final cluster centers.
  gamma=
    gammab[
      ,
      ((l1c-1)*KCn+1):
        (l1c*KCn)
    ]
  
  
  ## ==========================================================
  ## Recalculate beta and gamma using final clusters
  ## ==========================================================
  
  SS2=0
  
  beta=gamma=c()
  
  
  # Final observation-specific beta:
  # each observation is assigned the mean of its final cluster.
  
  for(i in 1:n)
    beta=c(
      beta,
      cols2(X[cxtn==cxtn[i],])/
        sum(cxtn==cxtn[i])
    )
  
  
  # Final cluster centers.
  for(i in 1:KCn)
    gamma=cbind(
      gamma,
      cols2(X[cxtn==i,])/
        sum(cxtn==i)
    )
  
  
  ## ==========================================================
  ## Calculate final within-cluster quadratic loss
  ## ==========================================================
  
  beta2=beta
  beta2=matrix(beta2,p,n)
  
  
  for(i in 1:KCn){
    
    # Singleton cluster.
    if(sum(cxtn==i)==1){
      
      SS2=
        SS2+
        sum(
          (
            sqrtm(solve(Csig))%*%
              (X[cxtn==i,]-gamma[,i])
          )^2
        )/(2*n)
      
    }else{
      
      # Cluster containing multiple observations.
      SS2=
        SS2+
        sum(
          (
            sqrtm(solve(Csig))%*%
              (t(X[cxtn==i,])-gamma[,i])
          )^2
        )/(2*n)
    }
  }
  
  
  ## ==========================================================
  ## Display result for the current lambda
  ## ==========================================================
  
  # Print:
  #   lambda  : current tuning parameter
  #   SS      : final within-cluster quadratic loss
  #   clusters: number of observations in each cluster
  
  cat(
    "lambda =", lambda,
    " SS =", SS2,
    " clusters =", paste(table(cxtn), collapse=","),
    "\n"
  )
  
  
  ## ==========================================================
  ## Store performance measures for the current lambda
  ## ==========================================================
  
  # Store the quadratic loss.
  S2=c(
    S2,
    sum(
      (
        sqrtm(solve(Csig))%*%
          (t(X)-beta2)
      )^2
    )/(2*n)
  )
  
  
  # Store mean squared error of beta.
  msel=c(
    msel,
    sum((beta-beta0)^2)/(n*p)
  )
  
  
  # Store cluster assignments.
  cxx=cbind(cxx,cxtn)
  
  
  # Store beta estimates.
  betat=cbind(betat,beta)
  
  
  # Store all gamma estimates generated during optimization.
  gammat=cbind(gammat,gammab)
}


## ============================================================
## Select the best lambda / final model
## ============================================================

# Select the model corresponding to the smallest objective
# function value stored in SS.

cxp=
  cxx[
    ,
    which.min(SS)
  ]


# Extract the corresponding beta estimate.
betap=
  betat[
    ,
    which.min(SS)
  ]


# Extract the corresponding cluster centers.
gammap=
  gammat[
    ,
    ((which.min(SS)-1)*KCn+1):
      (which.min(SS)*KCn)
  ]


## ============================================================
## Transform cluster centers back to original scale
## ============================================================

# X was standardized by:
#
#   X_original %*% solve(sqrtm(Csig2))
#
# Therefore, transform the estimated gamma back to the
# original data scale by multiplying sqrtm(Csig2).

gammap=
  sqrtm(Csig2)%*%
  gammap


## ============================================================
## Restore original data matrix
## ============================================================

X=
  Xlg[
    ((j-1)*n+1):(j*n),
    :
  ]


## ============================================================
## Store final results
## ============================================================

# Store the selected cluster assignment.
c_1210=
  rbind(
    c_1210,
    cxp
  )


# Store the selected cluster centers on the original
# data scale.
gamma_1210=
  rbind(
    gamma_1210,
    gammap
  )