Cessation of a salmon decline with control of parasites
================
Stephanie Peacock
2022-06-28

Stephanie J. Peacock, Martin Krkošek, Stan Proboszcz, Craig Orr, and
Mark A. Lewis. 2013. Cessation of a salmon decline with control of
parasites. Ecological Applications 23:606–620.
<http://dx.doi.org/10.1890/12-0519.1>

Details of the analysis of pink salmon population data, including
spawner and catch data and R code for compiling spawner-recruit pairs
and fitting the Ricker model.

## File list

`nuSEDS_PINK.csv`

`Catch_PINK.csv`

`1_DataCompilation.R`

`2_PopulationAnalysis.R`

## Description

The following files are made available with the intention of ensuring
our manipulation of data and analysis of pink salmon populations is
entirely transparent. We assume that the reader is familiar with the
program R, which can be downloaded from the web at cran.r-project.org.
Any errors in transcribing and manipulating the data are the sole
responsibility of the authors. Additional data and code for the analysis
of farm treatments and lice on wild salmon are available upon email
request to Stephanie Peacock (stephanie.j.peacock at gmail.com).

`nuSEDS_PINK.csv` contains spawner data for pink salmon populations in
British Columbia, Canada from 1950-2010, provided by Fisheries and
Oceans Canada (contact: Bruce Baxter <bruce.baxter@dfo-mpo.gc.ca>).
These data include the following columns:

-   Area: Fisheries management area, the spatial scale at which catch
    data are reported.
-   River: The river at which spawners were enumerated, the finest scale
    at which spawner data are available.
-   Yr: The year spawners were enumerated.
-   Escapement: The estimated number of spawners in that river in that
    year, from the nuSEDS database.

`Catch_PINK.csv` contains catch data for pink salmon by management area,
provided by different area managers at Fisheries and Oceans Canada
(contacts: Pieter VanWill <pieter.vanwill@dfo-mpo.gc.ca> and David
Peacock <david.peacock@dfo-mpo.gc.ca>). These data include the following
columns:

-   Area: Fisheries management area, the spatial scale at which catch
    data are reported.
-   Odd_Even: Integer indicated whether the catch was of an odd-year
    population (=1) or even year (=2).
-   Year: The year of the catch.
-   Catch: The number of pieces of pink salmon caught in all fisheries
    for the given year and area.

`1_DataCompilation.R` is R code that calls the previous two data files
and calculates the number of recruits per spawner.

`2_PopulationAnalysis.R` is R code that fits the linearized Ricker model
to the spawner recruit data, and tests for an effect of sea lice on wild
salmon.
