function [scenario_name, signals, init, sim_params] = scenario_2D_constant_wind_massless_controller()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~, signals, init, sim_params] = scenario_2C_constant_wind();

sim_params.use_massless_controller = true;
scenario_name = mfilename;

end