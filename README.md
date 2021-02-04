# WildBoarIPM
Code associated with formulating, implementing, and visualizing an analysis of wild boar population dynamics using an integrated population model (IPM).

## File description

 * `WildBoarIPM.R` <br/>
 Code for formulating, implementing, and running the model in NIMBLE.  <br/>
 Running requires an additional data file (`WildBoarIPM_Data.RData`). See below for information on data availability.
 
 <br/>
 
 * `WildBoarIPM_InitialValuesSim.R` <br/>
 Code for simulating initial values for running the IPM. 
 
 <br/>
   
 * `WildBoarIPM_Plots.R` <br/>
 Code for visualizing model results and generating the component plots for the figures included in the manuscript. <br/>
 Running requires posterior samples from the run model (`WildBoarIPM_MCMCsamples.RData`), which can be downloaded from: [TBA]
 <br/>
 
 ## Data availability
 
 To run the model, a file containing processed individual-based data on wild boar demography is required (`WildBoarIPM_Data.RData`). <br/>
 At present, this data is not openly available.  <br/>
 To request access to the data, please contact Dr. Marlène Gamelon: marlene.gamelon@univ-lyon1.fr
