# NestedPSO Algorithm for Finding Standardized Maximin *D*-optimal Designs for Enzyme Inhibition Models

This repository consists of the source codes of the online app, [Application of Particle Swarm Techniques in Finding Locally D-optimal and Standardized Maximin D-optimal Designs for Enzyme Inhibition Models](https://pingyangchen.shinyapps.io/stdmmoptdesigninhibition/).

We provide the source codes here for users who have basic knowledge in R programming.  There are two required R packages for running our code, *Rcpp* and *RcppArmadillo*, because the main function of NestedPSO is written in C++.  Nonetheless, users are NOT required to be familiar with C/C++ programming.  

In our online app, we search minimally supported standardized maximin *D*-optimal design only to avoid the huge computational time (hours to days, depending on the settings of swarm sizes and iterations of PSO loops) of one implementation that could exhaust the computational resource on the server.  Therefore, we suggest users to run the code in their own device if desired to change the number of support points.

The following R codes (*run.R*) are an example for how to use the source codes.  For more details of background knowledge, please refer to the main page of our online app or our published paper.

* Chen P.-Y., Chen, R.-B., Tung, H.-C., and Wong, W. K. (2017+). Standardized Maximim *D*-optimal Designs for Enzyme Kinetic Inhibition Models. (*Submitted*)

```R
### FILENAME: run.R
## Install Rcpp packages
install.packages(c("Rcpp", "RcppArmadillo"))
## Load Rcpp packages
library(Rcpp)
library(RcppArmadillo)

## Load Main Functions
sourceCpp("/kernel/psokernel.cpp")
source("rOptDesignInhibitionPSO.R")

## Get Target Design Information
D_INFO <- getD_INFO(model = 4, 
                    paraUpper = c(1, 5, 3, 5), 
                    paraLower = c(1, 4, 2, 4))
## Check the contents in D_INFO
names(D_INFO)
# If desired to change configuration, assign new values with '$' and variable lable
# For example, to change the upper bound of design space,
# D_INFO$dsUpper <- c(15, 90) # upper bound of c(s, i)
# To increase the number of support points,
# D_INFO$nSupp <- 5 

## Get PSO Configurations
PSO_INFO <- getPSO_INFO()
## Check the contents in PSO_INFO
str(PSO_INFO)
# If desired to change configuration, assign new values with '$' and variable lable
# For example, to increase the swarm size for OUTER loop,
# PSO_INFO$OUTER$nSwarm <- 200

## Run NestedPSO (PSO) Code
RESULT <- rOptDesignInhibitionPSO(D_INFO, PSO_INFO)
## Draw plot of directional derivative function
drawEquiv(RESULT, D_INFO)
```



