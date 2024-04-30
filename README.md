#wsyn: Wavelet approaches to studies of synchrony in ecology - introduction to the repository of codes supporting the talk

Daniel C. Reuman, University of Kansas, reuman@ku.edu  
Jonathan Walter, University of Virginia, jawalter@ucdavis.edu  

## Introduction

This repository can be used to reproduce the analyses behind the talk
``wsyn: Wavelet approaches to studies of synchrony in ecology'' given by 
Jon Walter and Dan Reuman as part of the ESA/Ecological Forecasting Initiative 
Webinar, 5/6/24, 12:00pm EDT. The idea behind this repo is to provide the 
complete, reproducible workflow for all the analysis and simulations that went 
into the talk, as well as the latex which produced the talk slides themselves.

## How to reproduce the workflow

To run the codes that reproduce the analyses of the paper, make your R working directory equal to the `code` directory of 
the repository 
and run `source("main.R")`. If all dependencies are in place (see below), this will pull data from the `data` directory and create 
results and then put them in the `Results` directory. 
You may have to install an R version or some R packages
that you don't have for this to work. Keep in mind these codes were run shortly prior to when the talk was given,
on Jon's and Dan's  machines, so no guarantees if a lot of time has passed since then.
We made no efforts at long-term reproducibility, but our R and package versions are listed below.

To recompile the talk slides, compile the latex in `talk`, under the file name `TheTalk.tex`.
This is ordinary beamer latex. 

## Dependencies

### Overview of dependencies

R, R studio, a setup that can compile beamer latex documents. The right versions of everything may be required - see below.

### Dependencies on R, R studio

Dan used R version 4.3.0 running on Ubuntu 18.04, and R studio version 2023.03.0+386. 
Jon used <Jon fill in your info here>

### Package versions

The packages we used were more or less current around the time the talk was given. Specifically:

Dan: 
wsyn, version 1.0.4;

Jon:
<Jon fill in yoru packages here>

You may have to install the versions of R, R studio, and the packages listed above to get the code to
run if a lot of time has passed since the talk. 

### Additional dependencies

The code and paper compilation process was tested by Jon and Dan on their machines. 
It has not been tested on other machines. We have endeavored to list all dependencies we can think of above, but we have 
only run the code on our machines, so we cannot guarantee that additional dependencies will not also be needed on other 
machines. This repository is intended to record a workflow, and is not designed or tested for distribution and wide use on 
multiple machines. It is not guaranteed to work on the first try without any hand-holding on arbitrary computing setups.

## Acknowlegements

<Insert ack here>













