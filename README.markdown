## Overview
 I've called the package 'spm' but the directory is named matrixmodels

Hopefully you can get the package installed using devtools by 
1. setting the working directory to Dropbox
2. install('matrixmodels')


## Development instructions 

Edit and add new functions in the files in the R/ directory. 
For conceptually separable sets of functions, add a new file. 

For best results use Hadley Wickham's `devtools`. Note you must also install `roxygen2`. 

`install.packages('devtools')`
`load(devtools)`

For common tasks, etc, see the [devtools github page](https://github.com/hadley/devtools)

As a quick start, modify away in the R/ directory and then run install('matrixmodels') or load_all('matrixmodels') again in R to have updates your interactive session.

Perhaps the first thing to do is edit the DESCRIPTION file to add yourself as author. 


## Additional documentation

### `pva.R`

The _x_`.iid` functions are all designed to compute based on random identically distributed draws _of complete matrices_. 

Yet to implement: 

* randomness based on _vital-rates_ See chapter 8 of Morris and Doak. 


