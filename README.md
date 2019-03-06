# DART-ID
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SlavovLab/DART-ID_2018/master)

Code for figures/data analysis of the DART-ID publication:

* [https://biorxiv.org/content/early/2018/08/23/399121](https://biorxiv.org/content/early/2018/08/23/399121)

-----------------

While the DART-ID software is written in Python, all downstream analyses, as well as the first version of DART, are written in R. Some analysis and software development was also done with iPython notebooks.

All figures can be generated with the scripts in ```Rscripts/fig_scripts```

## Data

All data referenced in this code is available on [MassIVE (ID: MSV000083149)](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=ed5a1ab37dc34985bbedbf3d9a945535). ```evidence_updated.txt``` files referenced are from the alignment results, located in ```MSV000083149/other/Alignments``` (ftp://massive.ucsd.edu/MSV000083149/other/Alignments/). R Data files can be found at ```MSV000083149/other/Rdata``` (ftp://massive.ucsd.edu/MSV000083149/other/Rdata/). 

## Help

All R code was run on the following R session:

```
R version 3.4.4 (2018-03-15)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.14.1

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RColorBrewer_1.1-2 ggridges_0.5.0     viridisLite_0.3.0  bindrcpp_0.2.2     pracma_2.1.5       forcats_0.3.0      stringr_1.3.1     
 [8] dplyr_0.7.6        purrr_0.2.5        readr_1.1.1        tidyr_0.8.1        tibble_1.4.2       ggplot2_3.0.0      tidyverse_1.2.1   

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.18     cellranger_1.1.0 pillar_1.3.0     compiler_3.4.4   plyr_1.8.4       bindr_0.1.1      tools_3.4.4      jsonlite_1.5    
 [9] lubridate_1.7.4  nlme_3.1-137     gtable_0.2.0     lattice_0.20-35  pkgconfig_2.0.2  rlang_0.2.2      cli_1.0.0        rstudioapi_0.7  
[17] yaml_2.2.0       haven_1.1.2      withr_2.1.2      xml2_1.2.0       httr_1.3.1       hms_0.4.2        grid_3.4.4       tidyselect_0.2.4
[25] glue_1.3.0       R6_2.2.2         fansi_0.3.0      readxl_1.1.0     modelr_0.1.2     magrittr_1.5     backports_1.1.2  scales_1.0.0    
[33] rvest_0.3.2      assertthat_0.2.0 colorspace_1.3-2 utf8_1.1.4       stringi_1.2.4    lazyeval_0.2.1   munsell_0.5.0    broom_0.5.0     
[41] crayon_1.3.4 
```

Python code was run on Python 3.6.6 on an Anaconda3 distribution, with the following (relevant) packages:

```
anaconda                  2018.12                  py27_0  
anaconda-client           1.7.2                    py27_0 
conda                     4.5.12                   py27_0  
conda-build               3.17.6                   py27_0
cython                    0.29.2           py27h0a44026_0
ipykernel                 4.10.0                   py27_0  
ipython                   5.8.0                    py27_0  
ipython_genutils          0.2.0            py27h8b9a179_0  
ipywidgets                7.4.2                    py27_0 
jupyter                   1.0.0                    py27_7  
jupyter_client            5.2.4                    py27_0  
jupyter_console           5.2.0                    py27_1  
jupyter_core              4.4.0                    py27_0  
jupyterlab                0.33.11                  py27_0  
jupyterlab_launcher       0.11.2           py27h28b3542_0
matplotlib                2.2.3            py27h54f8f79_0  
networkx                  2.2                      py27_1 
numpy                     1.15.4           py27hacdab7b_0  
numpy-base                1.15.4           py27h6575580_0 
pip                       18.1                     py27_0
pytest                    4.0.2                    py27_0
scikit-learn              0.20.1           py27h27c97d8_0  
scipy                     1.1.0            py27h1410ff5_2 
setuptools                40.6.3                   py27_0 
tk                        8.6.8                ha441bb4_0
yaml                      0.1.7                hc338f04_2
```
