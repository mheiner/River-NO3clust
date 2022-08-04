---
output:
  html_document: default
  pdf_document: default
---
# River-NO3clust

## Code to accompany "Model-based clustering of trends and cycles of nitrate concentrations in Rivers across France"

### by Matthew Heiner, Matthew J. Heaton, Philip White, Benjamin Abbott, Camille Minuado, and Remi Dupas
-----------
## Preprocessing

Unprocessed data are located in the `data` folder. Prior to analysis, it is necessary to 

1. Install all required R packages (loaded in the .R scripts) prior to processing data or fitting the model.
2. Run `0_formatData0_mergeNO3_stations.R`, and 
3. Run `0_formatData1_covariates.R`. 
 
The above steps will create `MergedRiverData.csv` in the `data` folder, which is used by all other routines.

## Fitting the model
The script `2_fit_model.R` fits the clustering/estimation model via MCMC. The first 91 lines include user-supplied prior and MCMC settings that can be changed. We recommend running the model with multithreading by sourcing shell script `runmod.sh`. Posterior samples are saved as a `.RData` file in the `results` folder. The script `2_continue_MCMC.R` can be used to extend existing MCMC chains.

## Analysis
All other scripts in the primary folder are intended to be run interactively.
The files `GoF.R` (goodness-of-fit), `posteriorAnalysis.R`, and `posteriorClusteringSALSO.R` load saved posterior samples and walk through posterior analysis.

----------
## Data sources
Nitrate time series: http://www.naiades.eaufrance.fr/france-entiere#/

Departments of France: https://www.data.gouv.fr/fr/datasets/contours-des-departements-francais-issus-d-openstreetmap/

Elevation: https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1/view

Catchment area: https://www.data.gouv.fr/fr/datasets/bd-alti-r-75-m-250-m-1-000-m/

Lithology: http://www.geocatalogue.fr/Detail.do?id=6388

Hydrology (IDPR): http://www.geocatalogue.fr/Detail.do?id=13039

Hydrology (runoff, topography): http://geowww.agrocampus-ouest.fr/geonetwork/apps/georchestra/?uuid=518b3e0a-ee55-40cb-a3ed-da00e60505aa

CORINE Land Cover: https://www.geoportail.gouv.fr/donnees/corine-land-cover-2012

Eco-regions: https://www.eea.europa.eu/data-and-maps/data/biogeographical-regions-europe-3

