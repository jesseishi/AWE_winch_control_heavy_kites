function [scenario_name, signals, init, sim_params] = scenario_2_use_massless()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~, signals, init, sim_params] = scenario_1B_change_in_wind_speed();

scenario_name = mfilename;
sim_params.use_massless_controller = true;

end