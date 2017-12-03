# ZT_vs_b

Plots ZT vs. b and locate measurement data on the plot.

This function accepts measured thermoelectric data and places it on the maximum ZT vs. b map. The load line for the measured data against maximum ZT vs. b is also shown. 

## Inputs

The program takes measured thermoelectric quantities being either
1. Seebeck, Resistivity, and ZT
2. Seebeck, Resistivity, and Total Thermal Conductivity
3. Seebeck, Resistivity, and Lattice Thermal Conductivity
   
Notice that ZT and Thermal conductivity cannot be given to the program at the same time. 
The program will use one and calculate the other.

For details, see "Input parameters from experimental measurements" section below

## Method synopsis

1. Using the given Seebeck coefficient, the Fermi level is determined from a single parabolic
   band model. 
2. Using the given resistivity and the determined Fermi level, the resistivity and electronic 
   thermal conductivity at different Fermi energies are calculated.
3. Using the given ZT or thermal conductivity, all the rest thermoelectric quantities of interest
   can be calculated.
4. Vary the lattice thermal conductivity and locate the optimal Fermi energy that maximize the ZT
   for each lattice thermal conductivity value. Covert this into a ZT|max vs. b curve.

## Outputs

Figure 1: ZT_{max} vs. b_l
Figure 2: ZT_{max} vs. b_total
Figure 3: Seebeck coefficient vs. b
Figure 4: Electrical resistivity vs. b
Figure 5: Thermal conductivity vs. b
Figure 6: Lorenz number vs. b

This script follows the general procedure published in
"Transport property analysis method for thermoelectric materials: material quality factor and the effective mass model"
SD Kang, GJ Snyder - arXiv preprint arXiv:1710.06896, 2017

This function is written by

Authors: Xufeng Wang, Evan Witkosk, and Mark Lundstrom

Contact: lundstro@ecn.purdue.edu

Copyright 2017 Purdue University
