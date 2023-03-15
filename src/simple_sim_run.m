%% Set the parameters.
clear
close all
addpath('src')
addpath('src/helper')
[kite, tether, winch, environment] = load_params_mat("my_MegAWES", "parameters");
