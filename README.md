# Larger but younger fish when growth compensates for higher mortality in warmed ecosystem

This repo contains R code for analyzing individual growth data and catch per unit effort data in the heated Biotest Lake and surrounding reference area, to investigate how warming has affected growth, size structure and mortality of perch, and how that affects population size structure.

We fit linear and non-linear hierarchical Bayesian models using the R-package [brms](https://github.com/paul-buerkner/brms) to analyse growth and mortality rates using regression models, and the R-package [sizeSpectra](https://github.com/andrew-edwards/sizeSpectra) for fitting size spectra to catch data using maximum likelihood mehthods, following [*Edwards* et al. 2017](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12641/full) and [*Edwards* et al. 2020](https://www.int-res.com/abstracts/meps/v636/p19-33/)

**Authors:** [Max Lindmark](https://maxlindmark.netlify.app/), Malin Karlsson, [Anna G?rdmark](https://internt.slu.se/en/cv-originals/anna-gardmark/)

## How to replicate our analyses and navigate this repo

`data`
Only processed data ready for analysis are uploaded here for reproducing the results, please consult the authors before using. The raw data are available from database KUL for some of the years: https://www.slu.se/institutioner/akvatiska-resurser/databaser/kul/ and the rest is hosted by SLU. Individual data were collated in Huss et al (2019).

`R`
Contains code for analysis and data processing

`figures`
Contains figures of results

`output`
Contains .rds objects of model outputs due to long compuation times. Currently empty because files are huge...


