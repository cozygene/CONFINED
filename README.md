# CONFINED - CCA ON Features for INter-dataset Effect Detection
### An ultra-fast implementation for canonical correlation analysis (CCA)
*CONFINED* was developed for the purpose of capturing replicable sources of biological variability in methylation data. These sources include, for example, age, sex, and cell-type composition. Importantly, the variation captured by *CONFINED* does not include any variability from technical or batch effects.

*CONFINED* is implemented in R and requires packages **Rcpp** and **RcppArmadillo**. If you are using a mac and having installation issues, try installing homebrew or xcode prior to installing **Rcpp** and **RcppArmadillo**. If RcppArmadillo is still having installation issues try running:
```
mkdir ~/.R
echo PKG_LIBS=\"\" > ~/.R/Makevars
```

After installing the necessary packages, please ensure that you have also downloaded the file **rcppCCA.cpp** and put it in the same directory as **CONFINED.R**. **rcppCCA.cpp** contains the functions for performing CCA. Our CCA algorithm is entirely based on that of R package **CCA** by Gonzalez and Dejean. We simply translated their code into C++ code using packages from **RcppArmadillo**, and all credit for the algorithm goes to them.

## Usage
As input, *CONFINED* takes mandatory arguments:
- Two matrices of size _m_ by _n1_ and _m_ by _n2_ (the number of rows is the same but not necessarily the same number of columns), where _m_ > both _n1_ and _n2_
-  _t_ the number of features (methylation sites) to use
-  _k_ the number of components to save

The following inputs are optional:
-  _saveOP_ boolean, save the output (_CONFINED_ components for each dataset and the ranked list of features) or not 
-  _outfile_ the prefix for saving the output files (default is "Xi_CONFINED_components.txt" and "CONFINED_ranked_features.txt")
-  _thresh_ The threshold for determining the rank of the low-rank approximation of the input matrices in the feature-selection step of the _CONFINED_ algorithm. The default is .95, and if there are no canonical variables with correlation > _thresh_, the rank is set to 1.

*CONFINED* returns a list containing two items:
-  _X1comps_ the _k_ components for dataset1 produced by _CONFINED_ using _t_ features
-  _X2comps_ the _k_ components for dataset2 produced by _CONFINED_ using _t_ features


We provide two subsets of whole-blood methylation datasets from Liu et al.[1] and Hannum et al.[2] for an example of CONFINED's usage. 

First load the **CONFINED.R** file into R:
```
source("CONFINED.R")
```
_Note_ it is critical that the file **rcppCCA.cpp** is in the same directory as the working environment. If you wish to override this setting, please change line 9 in the **CONFINED.R** file to point to the proper directory and file. 

Then, load the datasets:
```
dat1<-read.table("demo_data1.txt")
dat2<-read.table("demo_data2.txt")
```

Run _CONFINED_, saving the output with prefix "demo":
```
results<-CONFINED(X1=dat1, X2=dat2, t=3000, k=10, outfile="demo")
```

_CONFINED_ will then save the files


You can also use the components to predict various sources of biological variability. We provide two files of cell-type proportion estimates from the reference-based algorithm of Houseman et al.[3] as an example:
```
cellests1<-read.table(file = "demo_cellcomp1.txt")
sapply(1:10, function(i) sapply(1:dim(cellests1)[2], function(j) summary(lm(cellests1[,j] ~ results$X1_comps[,1:i]))$r.squared) )

cellests2<-read.table(file = "demo_cellcomp2.txt")
sapply(1:10, function(i) sapply(1:dim(cellests2)[2], function(j) summary(lm(cellests2[,j] ~ results$X2_comps[,1:i]))$r.squared) )
```
The call to sapply will return a matrix where each entry _ij_ corresponds to the R^2 value of predicting the *i*th cell-type's proportion using 1:_j_ components. Performance may not be optimal in the case of the demo as we have only provided about 5% of the entire datasets for the purpose of efficiency.


## An ultra-fast implementation for canonical correlation analysis (CCA)
If you only wish to use this software for performing quick canonical correlation analysis, simply download the **rcppCCA.cpp** file and ensure you have installed packages **Rcpp** and **RcppArmadillo**. Afterwards, CCA can be run with the following code:
```
library(Rcpp)
library(RcppArmadillo)
sourceCpp("rcppCCA.cpp")
#General function for solving CCA
#Takes two matrices, X,Y of size mXn_1 and mXn_2, where m>n_1,n_2
CCA = function (X, Y){
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
  return(list(A= A, B= B, U = U, V = V))
}

```
_A_ and _B_ are the loadings that will project _X_ and _Y_ into maximally correlated space. _U=XA_ and _V=YB_ are defined as the canonical variables of _X_ and _Y_, where the columns of _X_ and _Y_ have been centered in order to make the columns of U orthogonal to each other as well as the columns of V orthogonal to each other. 



[1]Yun  Liu,  Martin  J  Aryee,  Leonid  Padyukov,  M  Daniele  Fallin,  Espen  Hesselberg,  Arni Runarsson,  Lovisa  Reinius,  Nathalie  Acevedo,  Margaret  Taub,  Marcus  Ronninger,  Klementy  Shchetynsky,  Annika  Scheynius,  Juha  Kere,  Lars  Alfredsson,  Lars  Klareskog,
Tomas  J  Ekstrom,  and  Andrew  P  Feinberg.   Epigenome-wide association  data  implicate dna methylation as an intermediary of genetic risk in rheumatoid arthritis. *Nature Biotechnology*, 31:142 EP –, 01 2013.
[2]Gregory  Hannum,  Justin  Guinney,  Ling  Zhao,  Li  Zhang,  Guy  Hughes,  SriniVas  Sadda, Brandy Klotzle, Marina Bibikova, Jian-Bing Fan, Yuan Gao, Rob Deconde, Menzies Chen, Indika Rajapakse, Stephen Friend, Trey Ideker, and Kang Zhang. Genome-wide methylation profiles reveal quantitative views of human aging rates. *Molecular  cell*, 49(2):359–367, 01 2013.
[3] Eugene  Andres  Houseman,  William  P.  Accomando,  Devin  C.  Koestler,  Brock  C.  Christensen, Carmen J. Marsit, Heather H. Nelson, John K. Wiencke, and Karl T. Kelsey.  Dna methylation arrays as surrogate measures of cell mixture distribution. *BMC Bioinformatics*, 13(1):86, May 2012.



# CONFINED
