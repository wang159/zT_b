# zT_vs_b

Plots zT vs. b and locate measurement data on the plot.

This function accepts measured thermoelectric data and places it on the maximum zT vs. b map. The load line for the measured data against maximum zT vs. b is also shown. 

This script follows the general procedure published in

**["Transport property analysis method for thermoelectric materials: material quality factor and the effective mass model"](https://arxiv.org/abs/1710.06896)
SD Kang and GJ Snyder - arXiv preprint arXiv:1710.06896, 2017**

## Inputs

The program takes measured thermoelectric quantities being either
1. Seebeck, Resistivity, and zT
2. Seebeck, Resistivity, and Total Thermal Conductivity
3. Seebeck, Resistivity, and Lattice Thermal Conductivity
   
Notice that zT and Thermal conductivity cannot be given to the program at the same time. 
The program will use one and calculate the other.

![alt text](https://preview.ibb.co/g9gXcw/input_param.png)

## Method synopsis

1. Using the given Seebeck coefficient, the Fermi level is determined from a single parabolic
   band model. 
2. Using the given resistivity and the determined Fermi level, the resistivity and electronic 
   thermal conductivity at different Fermi energies are calculated.
3. Using the given zT or thermal conductivity, all the rest thermoelectric quantities of interest
   can be calculated.
4. Vary the lattice thermal conductivity and locate the optimal Fermi energy that maximize the zT
   for each lattice thermal conductivity value. Covert this into a zT|max vs. b curve.

## Outputs

Figure 1: zT_{max} vs. b_l

![alt text](https://image.ibb.co/hQHr4b/fig1.png)

Figure 2: zT_{max} vs. b_total

![alt text](https://image.ibb.co/fNjM4b/fig2.png)

Figure 3: Seebeck coefficient vs. b

![alt text](https://image.ibb.co/k6Tg4b/fig3.png)

Figure 4: Electrical resistivity vs. b

![alt text](https://image.ibb.co/d91OHw/fig4.png)

Figure 5: Thermal conductivity vs. b

![alt text](https://image.ibb.co/bx5iHw/fig5.png)

Figure 6: Lorenz number vs. b

![alt text](https://image.ibb.co/f3DEPb/fig6.png)

This function is written by

Authors: Xufeng Wang, Evan Witkosk, and Mark Lundstrom

Contact: lundstro@ecn.purdue.edu

Copyright 2017 Purdue University
