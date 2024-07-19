Simulation Studies:

This folder contains simulation studies for spike detection analysis using various firing rate functions.

Simulation Scripts:

1. static_sim.R: simulates spike trains under a static firing rate function.
2. dynamic_sim.R: simulates spike trains under a dynamic firing rate function.
3. flat_sim.R: simulates spike trains with a constant firing rate over time 
but varying across trials.

Source Code required: 

1. RcppFunctions.cpp: contains Rcpp functions to estimate spike locations using time-varying 
penalty through the vanilla dynamic programming algorithm. 
Utilize Rcpp drastically enhances the computational speed.

2. simulate_sptrain.R: contains functions to simulate spike trains using static or dynamic or constant firing rates.

- simulate.BP: uses a static firing rate; 
- simulate.BP2: uses a dynamic firing rate; 
- simulate.BP3: uses a constant firing rate varying across trials.

For further information or questions, please refer to the documentation within each script file.
