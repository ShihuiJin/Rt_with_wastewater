This folder contains all the codes required for the manuscript Combining wastewater surveillance and case data in estimating the time-varying effective reproduction number.
data.R is the file for cleaning the clincial and wastewater datasets. It also sources dist.R, which contains distributions utilized in this study. 
data.simulation.R is the file for simulating data for validation purposes. 
models.R involves all the model variants utilized for this study (parameter estimation with Rstan). 
subsample.R encompasses the model variants for various subsampling strategies proposed in the manuscript. 
analysis.R includes codes for extracting the Rstan outputs. 
