### This is a work in progress. Please feel free to message me with questions ###
# EDM with Age Structured Data

This repository contains the data and code associated with the manuscript in review "Age structure augments the predictive power of time series for fisheries and conservation" by Tara E. Dolan, Eric P. Palkovacs, Tanya L. Rogers, and Stephan B. Munch.


- `analysis_functions.R` contains functions used in scripts for this project. 
- `make_simulation_data.R` does 100 runs of simulations 1-3. 
    - Inputs: 
        - simulated_data/predmatrix_shortosc.csv (age specific predation of predators on prey)
        - simulated_data/preymatrix.csv (age specific predation of prey on predator eggs)
        - analysis_functions.R
    - Outputs: 
        - simulated_data/
             - Simulation1_data.csv = the resulting simulated data.
             - Simulation2_data.csv = the resulting simulated data.
             - Simulation3_data.csv = the resulting simulated data.
             - sim1_info.csv = metadata about simulation 1 runs, including recruitment noise, etc. 
             - sim2_info.csv = metadata about simulation 2 runs, including recruitment noise, etc. 
             - sim3_info.csv = metadata about simulation 3 runs, including recruitment noise, etc. 
        - sim_diagnostics/
             - sim1_post_burnin.png = plot of each age class over time (mean over 100 runs) after burn in period. 
             - sim1_pre_burnin.png = plot of each age class over time (mean over 100 runs) before burn in period.
             - sim1_run1.png = plot of each age class over time (first run only) after burn in. 
             - sim1_ycvsitself.png = plot of each year class vs. itself one time step previously
             - sim1_ycvsprevyr.png = plot of each year class vs. previous age class (same time step). 
             - sim2_post_burnin.png = plot of each age class over time (mean over 100 runs) after burn in period. 
             - sim2_pre_burnin.png = plot of each age class over time (mean over 100 runs) before burn in period.
             - sim2_run1.png = plot of each age class over time (first run only) after burn in. 
             - sim2_ycvsitself.png = plot of each year class vs. itself one time step previously
             - sim2_ycvsprevyr.png = plot of each year class vs. previous age class (same time step). 
             - sim3_post_burnin.png = plot of each age class over time (mean over 100 runs) after burn in period. 
             - sim3_pre_burnin.png = plot of each age class over time (mean over 100 runs) before burn in period.
             - sim3_run1.png = plot of each age class over time (first run only) after burn in. 
             - sim3_ycvsitself.png = plot of each year class vs. itself one time step previously
             - sim3_ycvsprevyr.png = plot of each year class vs. previous age class (same time step). 
             
Analyses             

- `modelcomparison.R` does the analysis that compares hierarchical, mixed age, and single age, and total abundance models for the simulated data
    - Inputs: simulated_data/
        - Simulation1_data.csv
        - Simulation2_data.csv
        - Simulation3_data.csv
    - Outputs: modelcomparison_outputs/
        - modelcomparison_simulation_out.csv
- `tslength_ageclass.R` add description
    - Inputs: simulated_data/
        - Simulation1_data.csv
        - Simulation2_data.csv
        - Simulation3_data.csv
    - Outputs: tslength_ageclass_outputs/
        -crossfits1.csv
        -crossfits2.csv
        -crossfits3.csv
- `pairwise_rho.R` obtains pairwise rho values for both pearson and dynamic correlation between age classes for simulated data 
    - Inputs:  simulated_data/
        - Simulation1_data.csv
        - Simulation2_data.csv
        - Simulation3_data.csv
    - Outputs:  pairwise_rho_outputs/
        -lagged_correlation_supermatrices.csv 
        -lagged_correlation_sim_stats.csv
        -dyrho_sim_stats.csv
        -dyrho_supermatrices.csv
- TotalAbundance_comp.R  Takes a subset of the simulated time series comprise of 30 years and 20 age classes and compares fitting a model trained on an index of total abundance (summed total of all age classes prior to model input) to a hierarchical model fit to all 20 age classes with the resulting predictions for each age class summed a posteriori. 
      - Inputs: simulated_data/
        - Simulation1_data.csv
        - Simulation2_data.csv
        - Simulation3_data.csv
      - Outputs: modelcomparison_outputs/thirtypoints.csv
        
- `empirical_analysis.Rmd` fits models to empirical data and visualizes model fits. does the analysis that compares hierarchical, mixed age, and single age, and total abundance models for the empirical data. Does the pairwise rho analysis for the empirical data and visualises those fits. 
    - Inputs:  empirical_data/
        - CTtrawlabundanceage19872017.csv Connecticut DEEP Long Island Sound trawl survey (LISTS) index
        - SB_dataming.csv Biomass and numbers at age for striped bass from CTDEEP
        - KRFC.csv escapement data and river returns for Klamath River Fall Chinnook Salmon
    - Outputs:  empirical_data/
        - SB_fitstatsE9.csv model fits for striped bass hierarchical models fit with E=9       
        - KRFC_fitstatsE6.csv model fits for salmon hierarchical models fit with E=6
    - Outputs: supplementary_figs/
        - r2ETSbiom.png grid of Es and taus, measured with r-squared
        - r2ETKRFC.png same as above, but for salmon
        - RMSE_ETSBbiom.png grid of Es and taus, measured with RMSE
        - RMSE_ETKRFC.png same as above, but for salmon
    - Outputs: figures/
        - SBfitoverdataE9allAgg.png hierarchical, total abundance and aggregate fits for striped bass models over data, E=9
        - KRFCfitoverdataE5all.png hierarchical, total abundance and aggregate fits for salmon models over data, E=6
        
       
       
        

Visualization

- `SimObsPredvis.R` plots observed and predicted values for a hierarchical GP fit to simulated datasets.
   - Inputs: simulated_data/
        - Simulation1_data.csv
        - Simulation2_data.csv
        - Simulation3_data.csv
   - Outputs: figures/
        - Sim1obspreds_age.png
        -
        
misc

- `simulation_sandbox.R` a place to play with the simulations. Stores the old parameter values for the April 2022 draft of the manuscript.
- Simulation fitstats: (`Sim1fitstats_100data.csv`, `Sim2fitstats_100data.csv`, `Sim3fitstats_100data.csv`)
- hundred_data_comp.R is a analysis that was dropped from the publication. Compares performance for models trained on 100 data points allocated one of five ways. 
- 
