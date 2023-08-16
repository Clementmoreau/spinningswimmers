# Description

MATLAB code accompanying the publications

> *Generalised Jeffery's equations for rapidly spinning particles. Part 1: Spheroids* by M. P. Dalwadi, C. Moreau, E. A. Gaffney, K. Ishimoto, and B. J. Walker

> *Generalised Jeffery's equations for rapidly spinning particles. Part 2: Helicoidal objects with chirality* by M. P. Dalwadi, C. Moreau, E. A. Gaffney, B. J. Walker, and K. Ishimoto.

# Use

Minimal code for exploring the dynamics is contained in the "main" folder.

## "main" folder

Running either of the scripts `main_spheroidal.m` and `main_helicoidal.m` will simulate and plot the motion of a spheroidal or helicoidal particle in shear flow. The setup can be customised by editing the top section of each of these files. 

Plotting can be customised by editing `plot_solution.m` and `plot_sphere.m`, the latter of which is found in the `helpers` directory.

## "figures_part_I" and "figures_part_II" folders

Each script generates one of the figures from the publications above. 

## "movies" folder 

Scripts with "data" in their filename generate at MATLAB data file containing the data with which the corresponding movie is generated. It is necessary to run them before building the movies with the "frame" scripts.

Scripts with "frame" in their name build the movie using the generated data.

