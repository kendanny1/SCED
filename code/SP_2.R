## ============================================================
## Load required packages and initialize parameters
## ============================================================

rm(list=ls())

library(openxlsx)
library(expm)


## ============================================================
## Global parameters
## ============================================================

v=6                    # Simulation / dataset index
KCnn=7                 # Target / reference number of clusters
n=250                  # Number of observations per dataset


## ============================================================
## Initialize result storage
## ============================================================

result2=c()            # General result storage
risum=c()              # Store Rand Index for each dataset
result1=c()            # General result storage
KCs1=c()               # Store selected KC1 for each dataset
cxse=c()               # Store selected cluster assignments


## ============================================================
## Function for calculating column sums
## ============================================================

# If x is already a p-dimensional vector, return x directly.
# Otherwise, calculate the column sums of x.
cols2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=colSums(x)
  }
  return(out)
}


## ============================================================
## Loop over all 500 datasets
## ============================================================

# j_0102 indexes the simulation dataset.
#
# For each dataset, results generated under different KC1
# settings are compared and the best clustering solution
# is selected according to the within-cluster loss.

j_0102=1
  
  
## ==========================================================
## Initialize storage for the current dataset
## ==========================================================

# Store the loss value obtained from each candidate KC1.
SSS_1128=c()

# Store the corresponding cluster assignments.
cxk=c()

# Store the corresponding KC1 values.
KCs=c()

# Store the clustering evaluation values.
ris=c()


## ==========================================================
## First pass: search candidate KC1 from 10 to KCnn
## ==========================================================

# Check whether the corresponding result file exists.
#
# This loop loads the available model result files and
# obtains the estimated cluster centers and assignments.

for(KC1 in 10:KCnn){
  
  if(
    file.exists(
      paste(
        "/NA3/kendanny1/sep/kotzp15/hs_",
        KC1,"_",
        KCnn,"_",
        v,"_",
        n,"_",
        j_0102,
        ".RData",
        sep=""
      )
    )==TRUE
  ){
    
    
    ## --------------------------------------------------------
    ## Load model result and original data
    ## --------------------------------------------------------
    
    load(
      paste(
        "/NA3/kendanny1/sep/kotzp15/hs_",
        KC1,"_",
        KCnn,"_",
        v,"_",
        n,"_",
        j_0102,
        ".RData",
        sep=""
      )
    )
    
    load(
      paste(
        "/NA3/kendanny1/sep/kotzp15/data_",
        v,
        ".RData",
        sep=""
      )
    )
    
    
    ## --------------------------------------------------------
    ## Extract the estimated gamma and cluster assignment
    ## for the current simulation dataset
    ## --------------------------------------------------------
    
    SS2=0
    
    
    # Extract the estimated cluster centers corresponding
    # to dataset j_0102.
    gammap=
      gamma_1210[
        ((j_0102-1)*p+1):
          (j_0102*p),
      ]
    
    
    # Extract the estimated cluster assignment.
    cxp=
      c_1210[j_0102,]
    
    
    # Extract the observations corresponding to dataset
    # j_0102.
    X=
      Xlg[
        ((j_0102-1)*n+1):
          (j_0102*n),
      ]
    
    
    ## ========================================================
    ## Recalculate cluster centers from original X
    ## ========================================================
    
    gamma=c()
    
    # Calculate the mean vector of each estimated cluster.
    for(i in 1:KCnn){
      
      gamma=
        cbind(
          gamma,
          cols2(X[cxp==i,])/
            sum(cxp==i)
        )
    }
    
    
    ## ========================================================
    ## Estimate within-cluster covariance matrix
    ## ========================================================
    
    Csig3=matrix(0,p,p)
    
    for(i in 1:KCnn){
      
      # Special treatment for a singleton cluster.
      if(length(X[cxp==i,])==p){
        
        Csig3=
          Csig3+
          outer(
            X[cxp==i,]-gamma[,i],
            X[cxp==i,]-gamma[,i],
            FUN="*"
          )
        
      }else{
        
        # Within-cluster scatter matrix for clusters
        # containing multiple observations.
        Csig3=
          Csig3+
          (
            t(X[cxp==i,])-gamma[,i]
          )%*%
          t(
            t(X[cxp==i,])-gamma[,i]
          )
      }
    }
    
    
    # Convert the scatter matrix into a covariance matrix.
    Csig3=Csig3/n
  }
}


## ==========================================================
## Second pass: evaluate candidate KC1 solutions
## ==========================================================

# Traverse KC1 in reverse order.
#
# For each available result file, calculate the final
# within-cluster quadratic loss using the original data
# and the estimated cluster assignment.

for(KC1 in KCnn:10){
  
  if(
    file.exists(
      paste(
        "/NA3/kendanny1/sep/kotzp15/hs_",
        KC1,"_",
        KCnn,"_",
        v,"_",
        n,"_",
        j_0102,
        ".RData",
        sep=""
      )
    )==TRUE
  ){
    
    
    ## --------------------------------------------------------
    ## Load the corresponding model result
    ## --------------------------------------------------------
    
    load(
      paste(
        "/NA3/kendanny1/sep/kotzp15/hs_",
        KC1,"_",
        KCnn,"_",
        v,"_",
        n,"_",
        j_0102,
        ".RData",
        sep=""
      )
    )
    
    
    # Load the original simulated data.
    load(
      paste(
        "/NA3/kendanny1/sep/kotzp15/data_",
        v,
        ".RData",
        sep=""
      )
    )
    
    
    ## --------------------------------------------------------
    ## Extract the result for dataset j_0102
    ## --------------------------------------------------------
    
    SS2=0
    
    
    # Estimated cluster centers for the current dataset.
    gammap=
      gamma_1210[
        ((j_0102-1)*p+1):
          (j_0102*p),
      ]
    
    
    # Estimated cluster assignment for the current dataset.
    cxp=
      c_1210[j_0102,]
    
    
    # Original observations for the current dataset.
    X=
      Xlg[
        ((j_0102-1)*n+1):
          (j_0102*n),
      ]
    
    
    ## ========================================================
    ## Recalculate cluster centers
    ## ========================================================
    
    gamma=c()
    
    for(i in 1:KCnn){
      
      gamma=
        cbind(
          gamma,
          cols2(X[cxp==i,])/
            sum(cxp==i)
        )
    }
    
    
    ## ========================================================
    ## Calculate within-cluster quadratic loss
    ## ========================================================
    
    # Calculate:
    #
    #   (1 / 2n) *
    #   sum_i || Sigma^{-1/2}(X_i - gamma_{c_i}) ||^2
    #
    # This loss is used to compare the clustering solutions
    # obtained under different KC1 values.
    
    for(i in 1:KCnn){
      
      # If cluster i contains one observation,
      # X[cxp==i,] is a vector rather than a matrix.
      if(sum(cxp==i)==1){
        
        SS2=
          SS2+
          sum(
            (
              sqrtm(solve(Csig3))%*%
                (
                  X[cxp==i,]-gamma[,i]
                )
            )^2
          )/(2*n)
        
      }else{
        
        # For a cluster containing multiple observations,
        # calculate the matrix-based quadratic loss.
        SS2=
          SS2+
          sum(
            (
              sqrtm(solve(Csig3))%*%
                (
                  t(X[cxp==i,])-gamma[,i]
                )
            )^2
          )/(2*n)
      }
    }
    
    
    ## --------------------------------------------------------
    ## Store the loss and corresponding model information
    ## --------------------------------------------------------
    
    # Store the loss for this KC1.
    SSS_1128=
      c(
        SSS_1128,
        SS2
      )
    
    
    # Store the corresponding cluster assignment.
    cxk=
      rbind(
        cxk,
        cxp
      )
    
    
    # Store the corresponding KC1 value.
    KCs=
      c(
        KCs,
        KC1
      )
    
    
    # Store the Rand Index associated with the selected
    # model from the original optimization.
    ris=
      c(
        ris,
        ril[which.min(SS)]
      )
}


## ==========================================================
## Select the best clustering solution
## ==========================================================

# Only perform model selection if at least one valid
# candidate result exists.

if(length(SSS_1128)>0){
  
  
  ## --------------------------------------------------------
  ## Select the candidate with minimum within-cluster loss
  ## --------------------------------------------------------
  
  # Identify the clustering solution with the smallest
  # recalculated quadratic loss.
  cxs=
    cxk[
      which.min(SSS_1128),
    ]
  
  
  # Store the selected cluster assignment.
  cxse=
    rbind(
      cxse,
      cxk[
        which.min(SSS_1128),
      ]
    )
  
  
  ## ========================================================
  ## Calculate Rand Index for the selected clustering
  ## ========================================================
  
  rri=0
  
  
  # Compare every pair of observations.
  #
  # The Rand Index measures whether pairs of observations
  # are classified consistently between the estimated
  # clustering and the true clustering.
  
  for(i in 1:(n-1)){
    
    for(k4 in (i+1):n){
      
      rri=
        rri+
        (
          
          # Case 1:
          # Both methods assign the pair to the same cluster.
          
          (
            (cxs[i]==cxs[k4])*
              (
                Clg[
                  n*(j_0102-1)+i
                ]==
                  Clg[
                    n*(j_0102-1)+k4
                  ]
              )
          )+
            
            
            # Case 2:
            # Both methods assign the pair to different clusters.
            
            (
              (cxs[i]!=cxs[k4])*
                (
                  Clg[
                    n*(j_0102-1)+i
                  ]!=
                    Clg[
                      n*(j_0102-1)+k4
                    ]
                )
            )
          
        )/(n*(n-1)/2)
    }
  }
  
  
  ## --------------------------------------------------------
  ## Store evaluation results
  ## --------------------------------------------------------
  
  # Store the Rand Index of the selected clustering.
  risum=
    c(
      risum,
      rri
    )
  
  
  # Store the KC1 value corresponding to the selected
  # clustering solution.
  KCs1=
    c(
      KCs1,
      KCs[
        which.min(SSS_1128)
      ]
    )
  }
}


## ============================================================
## Save final model-selection results
## ============================================================

# Save the complete workspace containing:
#
#   - cxse  : selected cluster assignments
#   - KCs1  : selected KC1 values
#   - risum : Rand Index values
#   - other intermediate results
#
# version=2 is used for compatibility with older R versions.

save.image(
  paste(
    "hsm_",
    KCnn,"_",
    v,"_",
    n,"_",
    j_0102,
    ".RData",
    sep=""
  ),
  version=2
)