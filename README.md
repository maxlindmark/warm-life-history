# Larger but younger fish when growth compensates for higher mortality in warmed ecosystem

[![DOI](https://zenodo.org/badge/298284214.svg)](https://zenodo.org/badge/latestdoi/298284214)

[Lindmark, M.](https://maxlindmark.github.io/), Karlsson, M., and
[Gårdmark, A](https://internt.slu.se/en/cv-originals/anna-gardmark/).
2023, January 9. Larger but younger fish when growth outpaces mortality
in heated ecosystem. bioRxiv.
<https://www.biorxiv.org/content/10.1101/2022.04.13.488128v3> (Accessed
18 January 2023).

This repo contains R code for analyzing individual growth data and catch
per unit effort data in the heated Biotest Lake and the surrounding
reference area, to investigate how warming has affected growth, size
structure and mortality of perch.

We fit linear and non-linear hierarchical Bayesian models using the
R-package [brms](https://github.com/paul-buerkner/brms) to analyse
growth and mortality rates, and the R-package
[sizeSpectra](https://github.com/andrew-edwards/sizeSpectra) for fitting
size spectra to catch data using maximum likelihood methods, following
[*Edwards* et al.
2017](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12641/full)
and [*Edwards* et al.
2020](https://www.int-res.com/abstracts/meps/v636/p19-33/).

## How to replicate our analyses and navigate this repo

From scratch: `data/raw` The raw data are available from the database
KUL for some of the years:
<https://www.slu.se/institutioner/akvatiska-resurser/databaser/kul/> and
the rest is hosted by SLU. Note that this data need some work before it
is usable. Individual data were collated in [Huss et al
(2019)](https://doi.org/10.1111/gcb.14637). If you have any questions
about these data (e.g., translating column names or how to use it),
please, feel free to reach out to us. We suggest, however, that you
start directly from the scripts in `R/analysis` and use the already
cleaned data.

`data/cleaned` Contains merged and cleaned data, see scripts
00_age_length_key.Rmd and 00_process_catch_data.Rmd in `R/clean_data`).

`data/for_fitting` Contains data ready to go directly into model fitting
scripts (`R/analysis`). These files are created in scripts in
`R/clean_data` starting with 01-04 (1 = von Bertalanffy, 2 = growth size
scaling, 3 = catch curves, 4 = size-spectra).

`R/analysis` Contains R code for fitting models (using data in
`data/for_fitting`), plotting and creating figures.

`figures` Contains all main and supporting figures (`figures/supp`)

You can also download this repo and view knitted html files of all .Rmd
scripts, or e.g., start directly from the scripts in `R/analysis` and
use the already cleaned data.

## Main packages and version used:

"other attached packages" (from `sessionInfo()`)

[1] modelr_0.1.8\
[2] patchwork_1.1.1\
[3] RColorBrewer_1.1-3\
[4] tidybayes_3.0.1\
[5] tidylog_1.0.2\
[6] bayesplot_1.7.2\
[7] viridis_0.5.1\
[8] viridisLite_0.4.1\
[9] nlstools_1.0-2\
[10] brms_2.17.0\
[11] Rcpp_1.0.8\
[12] forcats_0.5.1\
[13] stringr_1.5.0\
[14] dplyr_1.0.10\
[15] purrr_0.3.4\
[16] readr_2.1.1\
[17] tidyr_1.2.0\
[18] tibble_3.1.8\
[19] ggplot2_3.3.6\
[20] tidyverse_1.3.2
