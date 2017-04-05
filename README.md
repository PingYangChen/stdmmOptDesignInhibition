# Standardized Maximin *D*-optimal Designs for Enzyme Inhibition Models

We demonstrate how to use PSO to find standardized maximin *D*-optimal design for 4 types fo enzyme inhibition models.  
This app runs on RShiny, please check if the users had shiny package installed.  If not, please install shiny package in R.  
Also, the nestedPSO algorithm is written by C++ including armadillo linear algebra library.  
Thus, the Rcpp related packages are needed.
```R
install.packages(c("shiny", "Rcpp", "RcppArmadillo"))
```
To use this app, please open your R software and send the following commands.
```R
library(shiny)
runGitHub("stdmmOptDesignInhibition", "PingYangChen")
```
To use this app, please follow the instruction below.

1. This app will take a few seconds to initiate NestedPSO kernel via Rcpp.
2. Choose the inhibition model in **Select One**.
3. Specify the model parameters:
 * The first parameter, **Vmax** is linear to the *D*-criterion. One can verify this by observing results from different Vmax values.
 * When all values in **Lower** equal to that in **Upper** , the PSO search for locally *D*-optimal design will be performed.
 * When at leat one value in **Lower** is smaller than that in **Upper**, the NestedPSO search for standardized maximin *D*-optimal design will be performed.
4. Specify the design space for each component, **[s]** and **[i]**.
5. Specify the PSO parameters:
 * For finding locally *D*-optimal design, **swarm size** and **maximal number of iterations** for the PSO (**Local Loop**) precedure are needed.
 * For finding the standardized maximin *D*-optimal design, the NestedPSO consists of two nested PSO loop, **Outer** and **Inner** loops. For each loop, please specify **swarm size** and **maximal number of iterations**. Note that, by default seeting, the NestedPSO takes about 40 seconds and it should be enough to find the optimal design.
 * Some miscellaneous options are available for altering, such as the cognitive and the social parameters, the velocity clamping constant and the descending mode of inertia weight.
6. Click on **Execute** button to start the NestedPSO search.

