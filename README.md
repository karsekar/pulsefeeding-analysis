# pulsefeeding-analysis

## Introduction
The data and MATLAB and Python code used for analysis and figure generation of *E. coli* pulsing/starvation project is provided here. All the MATLAB code was written in version 2015b.

## The parts
This Github is divided into eight sections:
* Part 1 - Data and analysis of lag time versus feeding frequency - contains all of the OD data and code used for the wild-type pulsing experiments. Code generates the following plots:
  * All OD versus time plots with linear threshold fitting 
  * Lag time versus feedrate
  * Feedrate versus linear growth rate
  * Total glucose needed for proliferation versus feedrate
  * Summary table
* Part 2 - Flow Cytometry data and analysis. Code compares counts versus OD.
* Part 3 - Ion data data and analysis for varying feedrates. Code parses the ion data and generates plots where feedrate is varied (0.06, 0.12, 0.18 mmol/g/h) with the negative control in the background (f = 0 mmol/g/h)
* Part 4 - Ion data data and analysis for antibiotic additions. Code parses the ion data and generates plots where f = 0.18 mmol/g/h and antibiotics (Chloramphenicol, Rifamycin, and AZT) are added.
* Part 5 - Code to find intersecting proteins amongst those involved in division and that are degraded.
* Part 6 - Data and code for perturbing lag time (e.g. protease inhibitor, FtsZ titration)
* Part 7 - Data and code for Crp KO pulsing experiments.
* Part 8 - FtsZ modeling code. Code used to generate and test the FtsZ model (numerical and analytical solution).
* Common - contains common functions used throughout the different parts.
