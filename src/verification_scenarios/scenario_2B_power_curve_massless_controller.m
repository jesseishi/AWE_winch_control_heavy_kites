function [scenario_name, signals, init, sim_params] = scenario_2B_power_curve_massless_controller()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~, signals, init, sim_params] = scenario_2A_power_curve();

sim_params.use_massless_controller = true;
scenario_name = mfilename;

end