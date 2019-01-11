#Download and install Rcpp and Rcpp Armadillo
message("Loading CONFINED, ultra-fast CCA")
#library(Rcpp)
#library(RcppArmadillo)
message("Loaded Rcpp, compiling cpp code")

#This is a C++ file with the code written for the method
#It should be available where this file was downloaded
#' @useDynLib CONFINED, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
#Rcpp::sourceCpp("./src/rcppCCA.cpp")
message("Successfully compiled cpp code")

#' General function for quickly performing canonical correlation analysis (CCA)
#'
#' This function performs CCA on two matrices. As input, it takes
#' two matrices, \eqn{X} and \eqn{Y},  of size \eqn{m} by \eqn{n_1} and
#' \eqn{m} by \eqn{n_2} respectively, where \eqn{m}>\eqn{n_1},\eqn{n_2} (i.e., both
#' matrices have the same number of rows, but not necessarily the same
#' number of columns). This code is based on the algorithm by
#' Gonzalez and Dejean from the package '\code{CCA},' we just simply translated
#' it to C++ using functions from RcppArmadillo for speed. The canonical variables are returned
#' in decreasing order or correlation.
#' @param X \eqn{m} by \eqn{n_1} matrix
#' @param Y \eqn{m} by \eqn{n_2} matrix
#' @return A  -  the loadings for \eqn{X}
#' @return B  -  the loadings for \eqn{Y}
#' @return U  -  canonical variables of \eqn{X}, calculated by column centering \eqn{X} and projecting it on \eqn{A}
#' @return V  -  canonical variables of \eqn{Y}, calculated by column centering \eqn{Y} and projecting it on \eqn{B}
#' @return cors - the correlations of each corresponding pair of canonical variables e.g. cor_\eqn{i} = corr(U[, \eqn{i} ], V[, \eqn{i} ])
#' @export
CCA = function (X, Y){
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  # Check that the input's valid
  if(dim(X)[1] != dim(Y)[1]){
    stop("Unequal number of rows in the datasets. Ensure that the matrices are in correct orientation,
          e.g. the matrix is not transposed. Else, try taking the union of the rows.")
  }
  if(dim(X)[2] > dim(X)[1]){
    stop("The number of columns cannot exceed the number of rows (decomposition will fail).")
  }
  tmp<-doRCCA(X, Y)
  #loadings that project the data into maximally correlated space
  A<-tmp[[1]]
  B<-tmp[[2]]

  #there are only d canonical varibles
  d = min(c(dim(A)[2], dim(B)[2]))
  # s = min(dim(Y_dataframe)[1],dim(X_dataframe)[1]) - 1.
  A = A[,1:d]
  B = B[,1:d]

  #return the canonical variables (after doing some centering)
  U = getUV(X, A)
  V = getUV(Y, B)
  return(list(A= A, B= B, U = U, V = V, cors = tmp[[3]]))
}

#CONFINED's feature selection step
get_ranked_features = function (X_, Y_, thresh_=.95){
  #First get a low rank form of the canonical variables
  tmp = CCA(X_,Y_)
  ncol_=sum(sapply(1:10, function(i) cor(tmp$U[,i], tmp$V[,i]))>= thresh_)
  if(ncol_ < 1){
    message("Select lower correlation threshold, all canonical variables
            have correlation < thresh")
    message("Setting rank to 1")
    ncol_=1
  }

  #Now get a low-rank, correlated form of the input matrices
  Xtilde = t(centerCPP(X_)) %*% tmp$U[,1:ncol_] %*% t(tmp$U)[1:ncol_, ]
  Ytilde = t(centerCPP(Y_)) %*% tmp$V[,1:ncol_] %*% t(tmp$V)[1:ncol_, ]

  #Get the correlation of the low-rank and raw input-matrices' rows
  X_feat_cors  = sapply(1:dim(tmp$U)[1], function(i) cor(X_[i,], Xtilde[,i]))
  Y_feat_cors  = sapply(1:dim(tmp$V)[1], function(i) cor(Y_[i,], Ytilde[,i]))

  #Average them, then sort them
  feat_cors         = (X_feat_cors + Y_feat_cors)/2.
  feat_cors_sort    = sort(feat_cors, index.return=T, decreasing = T)
  UV_T_ranked_index = feat_cors_sort$ix
  #Return the most correlated indices

  return(UV_T_ranked_index )
}

#' Perform the CONFINED algorithm
#'
#' Generates components that capture sources of variability that are
#' shared between two datasets.
#' @param X1 Input matrices. \eqn{m} by \eqn{n_1} matrix. Must have the same number of rows as \code{X2}
#' @param X2 Input matrices. \eqn{m} by \eqn{n_2} matrix. Must have the same number of rows as \code{X1}
#' @param t Number of rows (e.g. methylation sites) to use
#' @param k Number of CONFINED components to produce
#' @param thresh Correlation threshold for selecting the number of canonical variables to use when generating the low-rank approximations in the feature selection step
#' @param outfile Prefix for saving the results
#' @param saveOP Boolean flag for saving the components/feature ranks to a txt file
#' @return X1comps  -  \emph{CONFINED} components for \code{X1}. These capture shared variability between \code{X1} and \code{X2}
#' @return X2comps  -  \emph{CONFINED} components for \code{X2}. These capture shared variability between \code{X1} and \code{X2}
#' @export
CONFINED<-function(X1, X2, t, k, thresh=.95, outfile="", saveOP=TRUE){

  X1<-as.matrix(X1)
  X2<-as.matrix(X2)

  #Rank the features
  message("Ranking features...")
  feats<-get_ranked_features(X1, X2, thresh)

  #Do CCA using the top t features
  tmp<-CCA(X1[feats[1:t], ], X2[feats[1:t], ])

  #Generate CONFINED components
  X1_ccs<-t(X1[feats[1:t], ])%*%tmp$U
  X2_ccs<-t(X2[feats[1:t], ])%*%tmp$V

  if(saveOP){
    message("Writing files...")
    if(outfile == ""){
      write.table(X1_ccs, quote=F, sep='\t', file = paste0("X1_CONFINED_components_t_",t, ".txt"),
                  row.names = F, col.names = F)
      write.table(X2_ccs, quote=F, sep='\t', file = paste0("X2_CONFINED_components_t_",t, ".txt"),
                  row.names = F, col.names = F)
      write.table(rownames(X1)[feats], quote=F, sep='\t', file = paste0("CONFINED_ranked_features.txt"),
                  row.names = F, col.names = F)
    }else{
      write.table(X1_ccs, quote=F, sep='\t', file = paste0(outfile, "_X1_CONFINED_components_t_",t, ".txt"),
                  row.names = F, col.names = F)
      write.table(X2_ccs, quote=F, sep='\t', file = paste0(outfile, "_X2_CONFINED_components_t_",t, ".txt"),
                  row.names = F, col.names = F)
      write.table(rownames(X1)[feats], quote=F, sep='\t', file = paste0(outfile, "_CONFINED_ranked_features.txt"),
                  row.names = F, col.names = F)
    }

  }
  message("Done!")
  return(list(X1_comps = X1_ccs, X2_comps = X2_ccs))
}

