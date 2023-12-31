
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BUSULFAN APP

This repository contains two folders:

- *sim*: materials to run a simulation study to quantify the interest of
  therapeutic drug monitoring to manage the inter-individual variability
  of exposure to busulfan.  
- *app*: a shiny application for model-informed precision dosing of
  busulfan in Thai pediatric patients.

Both rely on [mapbayr](https://github.com/FelicienLL/mapbayr), a package
dedicated to MAP-Bayesian estimation of PK parameters in R.

# Disclaimer

The application was developped as a proof-of-concept of a model-informed
precision dosing software dedicated to the calculation of doses of
busulfan in the Thai pediatric population. It is meant to be used for
research, educational and communication purposes. It comes with no
guarantee. It is not an electronic medical device. The authors of this
application are not responsible for the use of this application in a
clinical context for the treatment of actual patients.

# Installation

## Prerequisites

The following packages, and their dependencies, are required to make the
application work:

``` r
install.packages("mapbayr") # and mrgsolve dependency
install.packages("shiny")
install.packages("ggplot2")
install.packages("tidyr") # and tibble and dplyr dependencies
install.packages("lubridate")
install.packages("rhandsontable")
```

If you are using a Windows computer, you need to have Rtools installed
to compile C++. Visit [the CRAN
website](https://cran.r-project.org/bin/windows/Rtools/) for more
information.

## Run the application

You can open the application from R by calling:

``` r
shiny::runGitHub("FelicienLL/BUSULFANAPP", subdir = "app")
```

Alternatively, you can clone the repository or download its content in a
.zip file after clicking on the green “Code” icon on the bottom right of
this web page. Unzip it and simply launch the application in the script
called “app/app.R”.

![](capture.png)
