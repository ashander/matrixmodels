## Overview
 I've called the package 'tamp' but the directory is named matrixmodels

Can hopefully install package using devtools
1. library(devtools)
2. install_github(repo='matrixmodels',username='ashander')

## Development instructions 

Edit and add new functions in the files in the R/ directory. 
For conceptually separable sets of functions, add a new file. 

For best results use Hadley Wickham's `devtools`. Note you must also install `roxygen2`. 

`install.packages('devtools')`
`load(devtools)`

For common tasks, etc, see the [devtools github page](https://github.com/hadley/devtools)

As a quick start, modify away in the R/ directory and then run install('matrixmodels') or load_all('matrixmodels') again in R to have updates your interactive session.

## Additional documentation

The _x_`.iid` functions are all designed to compute based on random identically distributed draws _of complete matrices_. 




