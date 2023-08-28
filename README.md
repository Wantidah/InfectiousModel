# InfectiousModel
 Population dynamics and infectious disease modelling of wild bovids  \
 Explore the population impact of infectious diseases with different traits 

 There are two R scripts: 
 1) script_models_&_plotting_thisone.R includes:
 * All models used in the main text 
 * Non-infectious disease model
 * Six Bovine infectious disease models with five different model structures:
    + Anthrax, SI
    + Bovine tuberculosis, SEI
    + Hemorrhagic septicemia, SIRS
    + Lumpy skin disease, SEIRS
    + Foot and mouth disease, SEIRMS/E
    + Brucellosis, SEIRMS/E
2) script_population_dynamic_5sp.R includes:
* Non-infectious disease model for five bovid species (gaur, banteng, water buffalo, mainland serow and Chinese goral)
* Each species has different demographic parameters.

The result files from our models can be download here, includes: 
1) df_ndiff_allmodels.csv:
 + For generating summary statistics and box plots.
2) pca_table.xlsx:
 + For generating PCA biplot
