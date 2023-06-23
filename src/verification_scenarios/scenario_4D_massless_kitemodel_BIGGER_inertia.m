function [scenario_name, signals, init, sim_params, winch] = scenario_4D_massless_kitemodel_BIGGER_inertia(winch)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~, signals, init, sim_params] = scenario_4A_massless_kitemodel();

scenario_name = mfilename;
winch.J_kgm2 = 1000*winch.J_kgm2;

end