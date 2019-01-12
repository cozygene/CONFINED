# CONFINED - CCA ON Features for INter-dataset Effect Detection
### An ultra-fast implementation for canonical correlation analysis (CCA)
*CONFINED* was developed for the purpose of capturing replicable sources of biological variability in methylation data. These sources include, for example, age, sex, and cell-type composition. Importantly, the variation captured by *CONFINED* does not include any variability from technical or batch effects.

*CONFINED* is implemented in R and requires packages **Rcpp** and **RcppArmadillo**. If you are using a mac and having installation issues, try installing homebrew or xcode prior to installing **Rcpp** and **RcppArmadillo**. 

You can simply install *CONFINED* using **devtools**
```
devtools::install_github("cozygene/CONFINED")
```
Please see troubleshooting at the bottom for compilation issues.


## Usage
As input, *CONFINED* takes mandatory arguments:
- Two matrices of size _m_ by _n1_ and _m_ by _n2_ (the number of rows is the same but not necessarily the same number of columns), where _m_ > both _n1_ and _n2_
-  _t_ the number of features (methylation sites) to use
-  _k_ the number of components to save (can save up to min{_n1_, _n2_ components})

The following inputs are optional:
-  _saveOP_ boolean, save the output (_CONFINED_ components for each dataset and the ranked list of features) or not 
-  _outfile_ the prefix for saving the output files (default is "Xi_CONFINED_components.txt" and "CONFINED_ranked_features.txt")
-  _thresh_ The threshold for determining the rank of the low-rank approximation of the input matrices in the feature-selection step of the _CONFINED_ algorithm. The default is .95, and if there are no canonical variables with correlation > _thresh_, the rank is set to 1.

*CONFINED* returns a list containing two items:
-  _X1comps_ the _k_ components for dataset1 produced by _CONFINED_ using _t_ features
-  _X2comps_ the _k_ components for dataset2 produced by _CONFINED_ using _t_ features


We provide two subsets of whole-blood methylation datasets from Liu et al.[1] and Hannum et al.[2] for an example of CONFINED's usage. 

First load the **CONFINED** packages into R:
```
library(CONFINED)
```

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
-  *demo_CONFINED_ranked_features.txt* A file containing the list of features as sorted by *CONFINED*'s feature selection step
-  *demo_X2_CONFINED_components_t_3000.txt* A file containing an _n_ by _
-  demo_X1_CONFINED_components_t_3000.txt

You can also use the components to predict various sources of biological variability. We provide two files of cell-type proportion estimates from the reference-based algorithm of Houseman et al.[3] as an example:
```
cellests1<-read.table(file = "demo_cellcomp1.txt")
sapply(1:10, function(i) sapply(1:dim(cellests1)[2], function(j) summary(lm(cellests1[,j] ~ results$X1_comps[,1:i]))$r.squared) )

cellests2<-read.table(file = "demo_cellcomp2.txt")
sapply(1:10, function(i) sapply(1:dim(cellests2)[2], function(j) summary(lm(cellests2[,j] ~ results$X2_comps[,1:i]))$r.squared) )
```
The call to sapply will return a matrix where each entry _ij_ corresponds to the R^2 value of predicting the *i*th cell-type's proportion using 1:_j_ components. Performance may not be optimal in the case of the demo as we have only provided about 5% of the entire datasets for the purpose of efficiency.


## An ultra-fast implementation for canonical correlation analysis (CCA)
If you only wish to use this software for performing quick canonical correlation analysis, it can be accessed from the *CONFINED* package:
```
CONFINED::CCA(X,Y)
```
_A_ and _B_ are the loadings that will project input matrices _X_ and _Y_ into maximally correlated space. _U=XA_ and _V=YB_ are defined as the canonical variables of _X_ and _Y_, where the columns of _X_ and _Y_ have been centered in order to make the columns of U orthogonal to each other as well as the columns of V orthogonal to each other. 

Our CCA algorithm is entirely based on that of R package **CCA** by Gonzalez and Dejean. We simply translated their code into C++ code using packages from **RcppArmadillo**, and all credit for the algorithm goes to them.


## Troubleshooting
OSX (Mac) may have a problem where "math.h" is not found. This can usually be mitigated by running in the Terminal:
```
xcode-select --install
```

Another common problem on OSX (Mac) is an error concerning "lgfortran" or "quadmath." Below, we list steps suggested by [The Coatless Professor](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/). On that website, there are invaluable troubleshooting steps. Here, we will attempt to give the smallest number of required steps to take. Please visit the link for further details.
### R >= 3.5.x
<details><summary>Instructions</summary>
Copy and paste this into your Terminal window:

``` 
########### Xcode CLI

# Headless install of Xcode CLI
# Based on a script by Timothy Sutton, MIT licensed 2013 - 2014
# The code used is given at:
# https://github.com/timsutton/osx-vm-templates/blob/ce8df8a7468faa7c5312444ece1b977c1b2f77a4/scripts/xcode-cli-tools.sh#L8-L14

# Check if the Xcode CLI tool directory exists.
# See technical note: https://developer.apple.com/library/content/technotes/tn2339/_index.html#//apple_ref/doc/uid/DTS40014588-CH1-WHAT_IS_THE_COMMAND_LINE_TOOLS_PACKAGE_
# Note: This is not a rigorous check... So, if a user has deleted contents
# inside the folder but left the folder intact, then this will _not_ trigger
# an installation
if [ ! -d "/Library/Developer/CommandLineTools" ]; then

  # Create a temporary file for the header
  touch /tmp/.com.apple.dt.CommandLineTools.installondemand.in-progress

  # Figure out the correct Xcode CLI for the given mac OS
  PROD=$(sudo softwareupdate -l |
    grep "\*.*Command Line" |
    tail -n 1 | awk -F"*" '{print $2}' |
    sed -e 's/^ *//' |
    tr -d '\n')

  # Install Xcode CLI    
  sudo softwareupdate -i "$PROD" --verbose;

  rm -rf /tmp/.com.apple.dt.CommandLineTools.installondemand.in-progress
else
  echo "Xcode CLI is installed..."  
fi

########### clang6

# Download and Install the clang6 binary 
# Download ~440mb -> 2 gb installed
curl -O https://cran.r-project.org/bin/macosx/tools/clang-6.0.0.pkg
sudo installer -pkg clang-6.0.0.pkg -target /
```
Enter your password, then enter:
```
rm -rf clang-6.0.0.pkg

# Create an R environment file if it doesn't exist to store a modified path
# VARIABLE
if [ ! -e "~/.Renviron" ] ; then
   touch ~/.Renviron
fi

# Add the clang6 binary path to R's local paths
echo 'PATH="/usr/local/clang6/bin:${PATH}"' >> ~/.Renviron

########### gfortran

# Download and install the gfortran used in R 3.5.0
curl -O https://cloud.r-project.org/bin/macosx/tools/gfortran-6.1.pkg
sudo installer -pkg gfortran-6.1.pkg -target /
```
Enter your password once more (if prompted), and lastly:
```
rm -rf gfortran-6.1.pkg

# Establish a symlink of gfortran into /usr/local/bin
sudo ln -s /usr/local/gfortran/bin/gfortran /usr/local/bin/gfortran
```

If the above does not work and you've upgraded from R 3.0.0-3.3.3, try removing the old gfortran build, then reinstall the latest gfortran build:
```
# Download installer into working directory
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2

# Remove _files_ associated with the binary
for file in $(tar tfz gfortran-4.8.2-darwin13.tar.bz2); do
sudo rm -f /$file; 
done

# Remove empty _folders_ associated with the binary
for file in $(tar tfz gfortran-4.8.2-darwin13.tar.bz2); do 
sudo rmdir -p /$file; 
done

# Delete the installer
rm -rf gfortran-4.8.2-darwin13.tar.bz2

# Run the above step again
curl -O https://cloud.r-project.org/bin/macosx/tools/gfortran-6.1.pkg
sudo installer -pkg gfortran-6.1.pkg -target /
rm -rf gfortran-6.1.pkg

# Establish a symlink of gfortran into /usr/local/bin
sudo ln -s /usr/local/gfortran/bin/gfortran /usr/local/bin/gfortran
```
</p>
</details>

### 3.4.x
<details><summary>Instructions</summary>
The same link from the 3.5.x section will still be of help. You may try installing <a href="https://github.com/rmacoslib/r-macos-rtools/releases/download/v1.1.0/macos-rtools-1.1.0.pkg">these tools</a> from The coatless professor.
</p>
</details>

### 3.0.0-3.3.x
<details><summary>Instructions</summary>

Detailed instructions are provided by The coatless professor [here](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-os-x-before-r-3.4.0/).
Open the terminal and make sure xcode and gcc are installed:
```
xcode-select --install
```
Choose "Install" and verify that it was installed:
```
gcc --version
```
Now type:
```
cd /Applications/Utilities
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```
</p>
</details>




[1]Yun  Liu,  Martin  J  Aryee,  Leonid  Padyukov,  M  Daniele  Fallin,  Espen  Hesselberg,  Arni Runarsson,  Lovisa  Reinius,  Nathalie  Acevedo,  Margaret  Taub,  Marcus  Ronninger,  Klementy  Shchetynsky,  Annika  Scheynius,  Juha  Kere,  Lars  Alfredsson,  Lars  Klareskog,
Tomas  J  Ekstrom,  and  Andrew  P  Feinberg.   Epigenome-wide association  data  implicate dna methylation as an intermediary of genetic risk in rheumatoid arthritis. *Nature Biotechnology*, 31:142 EP –, 01 2013.
[2]Gregory  Hannum,  Justin  Guinney,  Ling  Zhao,  Li  Zhang,  Guy  Hughes,  SriniVas  Sadda, Brandy Klotzle, Marina Bibikova, Jian-Bing Fan, Yuan Gao, Rob Deconde, Menzies Chen, Indika Rajapakse, Stephen Friend, Trey Ideker, and Kang Zhang. Genome-wide methylation profiles reveal quantitative views of human aging rates. *Molecular  cell*, 49(2):359–367, 01 2013.
[3] Eugene  Andres  Houseman,  William  P.  Accomando,  Devin  C.  Koestler,  Brock  C.  Christensen, Carmen J. Marsit, Heather H. Nelson, John K. Wiencke, and Karl T. Kelsey.  Dna methylation arrays as surrogate measures of cell mixture distribution. *BMC Bioinformatics*, 13(1):86, May 2012.



