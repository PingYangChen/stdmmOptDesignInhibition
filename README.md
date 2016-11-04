# StdMxMnDesign_Inhibition

This app runs on RShiny, please check if the users had shiny package installed.  If not, please install shiny package in R.  
Also, the nestedPSO algorithm is written by C++ including armadillo linear algebra library.  
Thus, the Rcpp related packages are needed.
```R
install.packages(c("shiny", "Rcpp", "RcppArmadillo", "inline"))
```
To use this app, please open your R software and send the following commands.
```R
library(shiny)
runGitHub("StdMxMnDesign_Inhibition", "PingYangChen")
```
Note that, it is normal that this app might take some time to initiate the C++ kernel functions.  Once the app is ready, please follow the instruction below:

1. Choose the interested inhibition model in **"Select One"**.
2. Specify the model parameters. The users can change the value of the first parameter, *a*, to verify that it is linear to the *D*-criterion. For the lower and upper bounds for the remaining parameters, please find the illustation below.
  * When some values in **"Lower"** is smaller than that in **"Upper"**, the app will search the standardized maximin *D*-optimal design among the specified parameter space.
  * When all values in **"Lower"** equal to that in **"Upper"**, the app will search the locally *D*-optimal design at the specified parameter point.
3. Specify the design space for each component *s* and *i*.
4. Set the PSO options:
  * For searching the standardized maximin *D*-optimal design, there are two PSO procedures.  First, specify the wanted swarm size and maximal iteration number for the **"Outer"** and **"Inner"** loops of design search procedure, the NestedPSO.  Then, one also need to specify these two settings for the "Mu Loop", which is another PSO procedure for obtaining the assistant design used in the equivalence theorem.  **Note that, the NestedPSO takes about 40 seconds if the user runs it with default settings, and, it should be enough to find the optimal design.**
  * For searching the locally *D*-optimal design, the swarm size and maximal iteration number for the PSO (**"Local Loop"**) precedure is needed.
  * Some miscellaneous options are available for altering, such as the cognitive and the social parameters, the velocity clamping constant and the descending mode of inertia weight
5. As the settings are done, press **"Execute"** to start the PSO search.

