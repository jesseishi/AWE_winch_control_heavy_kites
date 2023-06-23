function [scenario_name, signals, init, sim_params, winch] = scenario_4C_massless_kitemodel_smaller_radius(winch)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~, signals, init, sim_params] = scenario_4A_massless_kitemodel();

scenario_name = mfilename;
winch.r_m = 0.08*winch.r_m;

end