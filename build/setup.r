#
if(!("RcppArmadillo" %in% installed.packages())){
install.packages("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.6.100.0.0.tar.gz", repos=NULL, type="source")
}

if(!("devtools" %in% installed.packages())){
install.packages("devtools")
}

if(!("igraph" %in% installed.packages())){
library(devtools)
install_github("igraph/rigraph")
}

cranPackages = c(
# to connect to jupyter
"repr",
"IRdisplay",
"evaluate",
"crayon",
"pbdZMQ",
"devtools",
"uuid",
"digest",

# for actual work
"tidyverse",
"gridExtra",
"vegan",
"imputeMissings",
"compositions",
"tidyverse",
"broom",
"knitr",
"kableExtra",
"IRdisplay",
"MiRKAT",
"purrrlyr")

new.packages <- cranPackages[!(cranPackages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# bioconductor packages

source('http://bioconductor.org/biocLite.R')
if(!("rhdf5" %in% installed.packages())){
biocLite("rhdf5")
}

if(!("biom" %in% installed.packages())){
library(devtools)
install_github("joey711/biom")
}

if(!("dada2" %in% installed.packages())){
library(devtools)
install_github("benjjneb/dada2")
}

bioconductorPackages <- c(
"phyloseq",
"dada2",
"qvalue")



new.packages.b <- bioconductorPackages[!(bioconductorPackages %in% installed.packages()[,"Package"])]
if(length(new.packages.b)) biocLite(new.packages.b)
