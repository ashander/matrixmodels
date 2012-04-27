## Overview
 I've called the package 'tamp' but the directory is named matrixmodels

On Mac and Unix, can install package using devtools

````R
install.packages('devtools')
library(devtools)
install_github(repo='matrixmodels',username='ashander')
````

For windows, can download a zip of the package using button above. 

## Development instructions 

Edit and add new functions in the files in the R/ directory. 
For conceptually separable sets of functions, add a new file. 

On Mac/Unix for best results use Hadley Wickham's `devtools`. Note you must also install `roxygen2`. 
For common tasks, etc, see the [devtools github page](https://github.com/hadley/devtools)
As a quick start, modify away in the R/ directory and then run install('matrixmodels') or 
load_all('matrixmodels') again in R to have updates your interactive session.

On other platforms, develop and/or modify code and document using the ROxygen style. 
For ROxygen style, just follow what is currently in the R/*.R files.

## Additional documentation

The _x_`.iid` functions are all designed to compute based on random identically distributed draws _of complete matrices_. 




