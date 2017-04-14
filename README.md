# Standardized Maximin *D*-optimal Designs for Enzyme Inhibition Models

We demonstrate how to use PSO to find standardized maximin *D*-optimal design for 4 types fo enzyme inhibition models.  
This app runs on RShiny, please check if the users had shiny package installed.  If not, please install shiny package in R.  
Also, the nestedPSO algorithm is written by C++ including armadillo linear algebra library.  
Thus, the Rcpp related packages are needed.

```R
# Install Rcpp packages
install.packages(c("Rcpp", "RcppArmadillo"))
# Load Rcpp packages
library(Rcpp)
library(RcppArmadillo)

path <- "D:/stdmmOptDesignInhibition"

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



