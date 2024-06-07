
static_sim.R/dynamic_sim.R generate results for simulated data in Tong et al. (2021)
https://arxiv.org/abs/2103.03818
The code are written in the way to run in server.


source code needed: function_l0spike.R; simulate_sptrain.R

function_l0spike.R: functions needed to estimate spikes using time-varying penalty
(main function: estspike.gaussian)

simulate_sptrain.R: functions needed to simulate spike trains using static/dynamic 
firing rates in Tong et al. (2021)
(simulate.BP: using static firing rate; simulate.BP2: using dynamic firing rate)