# EDM with Age Structured Data

This repository contains the data and code associated with the paper "Age structure augments the predictive power of time series for fisheries and conservation" by Tara E. Dolan, Eric P. Palkovacs, Tanya L. Rogers, and Stephan B. Munch.

**Tara, I recommend organizing the files into subfolders, so it's easier to find things. You will need to update the input/output file paths in the scripts (I did some but not all). Delete anything that is old or otherwise needs to be deleted (e.g. no longer in paper, file vs fileNEW, etc.). For stuff you want to keep around but not put on github, make a subfolder called 'archive' and put the files in there. I have listed 'archive' in the gitignore, so the files will be deleted from the online repo, but will still be retained locally on your computer. If you have questions, let me know.**

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
        - figures/
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
- `tslength_ageclass.R` add description (**I think this is parallelCrossfit.R? Please rename as appropriate.**)
    - Inputs:
    - Outputs:
- `pairwise_rho.R` obtains pairwise rho values between age classes for simulated and empirical data (**I think this is dyrho_parallel.R? Please rename as appropriate.**)
    - Inputs:  
    - Outputs:  
- `empirical_analysis.Rmd` does the analysis that compares hierarchical, mixed age, and single age, and total abundance models for the empirical data
    - Inputs:  
    - Outputs:  

Visualization

- `SimObsPredvis.R` plots observed and predicted values for a hierarchical GP fit to simulated datasets.
   - Inputs: simulated_data/
        - Simulation1_data.csv
        - Simulation2_data.csv
        - Simulation3_data.csv
   - Outputs: figures/
        - ?
- etc.

## Previous structure (delete)

- dyrho_parallel.R calculates the lagged correlation and pairwise dynamic rho for each simulation. 
=======
##### Updates 9/6/22 here ###
- simulation_sandbox.R = a place to play with the simulations. Stores the old parameter values for the April 2022 draft of the manuscript and also provides space for new values to try. 
    - Inputs: analysis_functions.R
- analysis_functions.R = a source for functions used in this project
- make_simulation_data.R = runs each simulation and produces the preylist files. 
    - Inputs: analysis_functions.R
    - Outputs: 
        - Simulation1_data.csv - data frame of 100 simulation runs of simulation 1
        - sim1_info.csv
        - Simulation2_data.csv
        - sim2_info.csv
        - Simultion2_data.csv
        - sim3_info.csv
- mixedAgeloop.R = compares single age, total abundance and hierarchical aged models
        - Inputs: Simulation1_data.csv, Simulation2_data.csv, Simulation3_data.csv
        - Outputs: mixedAgeOutnew.csv
- SimObsPredVis.R = makes the plots of data from 1st (of 100) run of the simulation with model over it for the simulations. 
        - Inputs: Simulation1_data.csv, Simulation2_data.csv, Simulation3_data.csv
        - Outputs:
            - Sim1obspreds_age.png a graph of the hierarchical model overlaying 30 years of data from one simulation run. (For the appendix). 
            - Sim1obspreds_agg.png  a graph of the hierarchical and total abundance model overlaying 30 years of data from one simulation run. 
            - Sim2obspreds_age.png a graph of the hierarchical model overlaying 30 years of data from one simulation run. (For the appendix). 
            - Sim2obspreds_agg.png  a graph of the hierarchical and total abundance model overlaying 30 years of data from one simulation run. 
            - Sim3obspreds_age.png a graph of the hierarchical model overlaying 30 years of data from one simulation run. (For the appendix). 
            - Sim3obspreds_agg.png  a graph of the hierarchical and total abundance model overlaying 30 years of data from one simulation run. 
            - Aggregatepreds20ages.png is a graph of the total abundance data, hierarchical model fit and total abundance model fit from one simulation run. 
- hundred_data_comp.R = splits the data into 5 ages, 20 years, 4 ages 25 years, 10 ages 10 years and 20 ages 5 years. Hierarchical model of each. Total abundance model of each, single age model of each.
        - Inputs: analysis_functions.R, Simulation1_data.csv, Simulation2_data.csv, Simulation3_data.csv
        - Outputs: Sim1Phis_100data.csv, Sim2Phis_100data.csv, Sim2Phis_100data.csv
        - Outputs: Sim1fitstats_100data.csv, Sim2fitstats_100data.csv, Sim3fitstats_100data.csv
- thirty_ptscomp.R = does the 30 years and 15 years of 20 age classes from all simulations total abundance, overall hierarchical fit and summed hierarchical fit comparison. Ultimately we did not include this analysis as we felt it was redundant to the hundred data comparison. 
        - Inputs: Simulation1_data.csv, Simulation2_data.csv, Simulation3_data.csv
        - Outputs: thirtypoints.csv, fifteenpoints.csv
- Sim_Data_vis100.Rmd = creates the visualizations for the 100 data points comparisons and the figures tanya requested.
    - Inputs: Sim1Phis_100data.csv, Sim2Phis_100data.csv, Sim2Phis_100data.csv
    - Outputs: phiplot.png - a graph of the mean value of phi from the hierarchical model from the hundred data comp script. 
    - Inputs: thirtypoints.csv or fifteenpoints.csv
    - Outputs: “thirtypointsfig.png”, thirty points comparison figure. 30 points of 20 ages, all three simulations r2 or RMSE. There are also the “thirtypointsI.png” “thirtypointsII.png” and “thirtypointsIII.png”, which are separated out for the different simulations. Ultimately decided not to include this figure. 
    - Inputs: im1fitstats_100data.csv, Sim2fitstats_100data.csv, Sim3fitstats_100data.csv
    - Outputs: 
        - cplotattempt4sumhier.png - a grob plot of the summed hierarchical fits compared to total abundance from 100 data comparison
        - “Boxplot_Hier_TA.png” - a box plot of the overall hierarchical fits compared to total abundance from 100 data comparison
        - “Boxplot_sumhier_TA.png” - a box plot of the individual hierarchical fits summed to total abundance compared to total abundance from 100 data comparison
        - “Boxplot_hier.png” - a box plot of the overall hierarchical fits vs time series length. 
- empirical_analysis.Rmd = does the empirical data analysis, the dynamic correlation plots for the empirical data and some visualizations of empirical data only. Outputs the csv that is used for the mixed analysis visualization. 

    - Inputs: 
        - "CTtrawlabundanceage19872017.csv" = striped bass abundance index
        - "SB_dataming.csv" = biomass and numbers at age
    - Outputs: 
        - "SB_rawData.png"  - time series of Striped bass by age
        - "r2ETSBbiom.png” & "RMSE_ETSBbiom.png”= Tile plots of SB fits across a grid of E and Tau using R-squared and RMSE respectively.
        - Lagged correlation and dynamic correlation plots - does not produce a .png file, but the plot may be saved manually. 
        - "SBfitoverdataE9allAgg.png" = Hierarchical fits over data 
        - “SB_fitstatsE9.csv” = total abundance, individual age and hierarchical fits 
    - Inputs: 
        - “KRFC.csv” = abundance at age of Klamath River fall chinook salmon. 
    - Outputs: 
        - "KRFC_rawData.png" = time series of KRFC abundance by age
        - Lagged correlation and dynamic correlation plots - does not produce a .png file, but the plot may be saved manually. 
        - “r2ETKRFC.png" & "RMSE_ETKRFC.png" - lagged fits across a grid of E and Tau using R-squared and RMSE respectively.
        - “KRFC_fitstatsE6.csv” = total abundance, individual age and hierarchical fits 
        - "KRFCfitoverdataE6all.png" = hierarchical fits over data
    - Outputs (combined): 
        - “mixedage_fitstats_emp.csv” and “mixedage_fitstats_empLOG.csv”  = combined csv of fits for data visualization, unlogged and logged version, respectively.  
- mixedAgeVis.Rmd - visualizations for the mixed age analysis
    - Inputs: "mixedage_fitstats_empLOG.csv”; "mixedageoutnew.csv"
    - Outputs: 
        - "MixedAge.png" - mixed age plot in 5 panels with vertical columns of points
        - “MixedAge1.png”, MixedAge2.png, MixedAge3.png and MixedAgeBottom.png are individual panels of the MixedAge.png plot that can be assembled in InDesign. 
        - Tinyviolins plots. The overall plot and the individual panels which can be assembled separately. 
- parallelcrossfit.R = performs the analysis which fits a hierarchical model to every possible combination of 1-20 age classes and 15 to 30 years. This analysis is best run on a multicore server because it can take a while, but it can be run on a standard machine. 
    - Inputs: Simulation1_data.csv, Simulation2_data.csv, Simulation3_data.csv
    - Outputs: crossfits_sim1.csv, crossfits_sim2.csv, crossfits_sim3.csv
- crossfitsViz.Rmd - visualizations for the parallel crossfits analysis
    - Inputs: crossfits_sim1.csv, crossfits_sim2.csv, crossfits_sim3.csv
    - Outputs: the parallel crossfits contour plots. These cannot be saved using ggsave. Usually I save them via screenshot.
- dyrho_parallel.R - does the comparison of dynamic correlation and lagged correlation for the simulated data. This analysis is best run on a multicore server because it can take a while, but it can be run on a standard machine. 
    - Inputs: Simulation1_data.csv, Simulation2_data.csv, Simulation3_data.csv
    -  
