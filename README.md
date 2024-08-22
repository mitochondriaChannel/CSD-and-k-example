# Cumulative Square Displacement (CSD)
The Cumulative Square Displacement is a novel quantitative tool that enables the detection of localized mechanical impulses acting on filamentous organelles that occur amidst the
stationary fluctuations caused by thermal jittering. We use this method to explore the forces affecting mitochondria within living cells; however, this approach can be expanded to other organelles or
biopolymers. The details of this can be found in "Deciphering the intracellular forces shaping mitochondrial motion" by Agustina Belén Fernández Casafuz, Azul María Brigante, María Cecilia De Rossi, Alejandro Gabriel Monastra and Luciana Bruno (currently in peer review). 

This repository contains the executable code necessary to obtain the CSD graph and calculate the number of events of active forces. The example corresponds to the time series of a mitochondrion in a *Xenopus laevis* melanocyte, shown in Fig. 3(b). 

The main script analyzes the data using the functions defined in **methods.py**. 

The input is the data file **snake_fig3b.txt**, which contains the coordinates of the mitochondrion shown in Fig. 3(b) for all the frames of the time-lapse. These coordinates were obtained by tracking the mitochondrion in confocal microscopy images using the ImageJ plugin JFilament (https://www.lehigh.edu/~div206/jfilament/). The raw file obtained from JFilament was preprocessed to convert the coordinates from pixels to $\mu$m and to have a time column as well. Then, the snake file displays 4 columns: frame number, time, X and Y coordinates. 

# License
Cumulative Square Displacements (CSD) © 2024 by Agustina Fernandez Casafuz, Azul Brigante, M. Cecilia De Rossi, Alejandro Monastra, and Luciana Bruno is licensed under CC BY-NC-SA 4.0
