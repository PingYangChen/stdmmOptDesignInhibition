# NestedPSO Algorithm for Finding Standardized Maximin *D*-optimal Designs for Enzyme Inhibition Models

This repository consists of the source codes of the online app, [Application of Particle Swarm Techniques in Finding Locally D-optimal and Standardized Maximin D-optimal Designs for Enzyme Inhibition Models](https://pingyangchen.shinyapps.io/stdmmoptdesigninhibition/).

We provide the source codes here for users who have basic knowledge in R programming.  There are two required R packages for running our code, *Rcpp* and *RcppArmadillo*, because the main function of NestedPSO is written in C++.  Nonetheless, users are not required to be familiar with C/C++ programming.  

By using the source codes, one can generate a design with more than minimally supported number of support points freely.  Once the standardized maximin design is not minimally supported, due to the lack of analytical solution of locally *D*-optimal design, the NestedPSO becomes a 3-layer optimization tool and the computational time would be exponentially increased.  Therefore, we 


```R
# Install Rcpp packages
install.packages(c("Rcpp", "RcppArmadillo"))
# Load Rcpp packages
library(Rcpp)
library(RcppArmadillo)

# Load Main Functions
sourceCpp("/kernel/psokernel.cpp")
source("rOptDesignInhibitionPSO.R")

# Get Target Design Information
D_INFO <- getD_INFO(model = 4, 
                    paraUpper = c(1, 5, 3, 5), 
                    paraLower = c(1, 4, 2, 4))
# Check the contents in D_INFO
names(D_INFO)
# If desired to change configuration, assign new values with '$' and variable lable
# For example, to change the upper bound of design space,
# D_INFO$dsUpper <- c(15, 90) # upper bound of c(s, i)

# Get PSO Configurations
PSO_INFO <- getPSO_INFO()
# Check the contents in PSO_INFO
str(PSO_INFO)
# If desired to change configuration, assign new values with '$' and variable lable
# For example, to increase the swarm size for OUTER loop,
# PSO_INFO$OUTER$nSwarm <- 200

# Run NestedPSO (PSO) Code
RESULT <- rOptDesignInhibitionPSO(D_INFO, PSO_INFO)
# Draw plot of directional derivative function
drawEquiv(RESULT, D_INFO)
```



