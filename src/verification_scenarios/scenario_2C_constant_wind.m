function [scenario_name, signals, init, sim_params] = scenario_2C_constant_wind()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
scenario_name = mfilename;

% One loop = 20 sec.
t_end = 40;
sim_params = struct('t_end', t_end, ...
    'use_massless_controller', false, ...
    'use_larger_winch_PI', false, ...
    'use_zgraggen', false, ...
    'use_massless_kitemodel', false);

% Set the windspeed.
t = (0:0.1:t_end)';
vw_mps = 10 * ones(size(t));
vw_mps(t>5) = 15;
vw_mps(t>10) = 20;
vw_mps(t>15) = 25;
vw_mps(t>20) = 30;

% vw_mps = linspace(10, 30, length(t))';
vw_mps = timeseries(vw_mps, t, 'Name', 'wind speed [m/s]');

% Other variables.
% These variables are variable in other verification scenario's so they
% must be compatible with `from_workspace`, meaning we must do this awkward
% construction of [constant, constant].
beta_deg = [30, 30];
phi_deg = [17.5, 17.5];
chi_unwrapped_deg = [90 90];
chi_deg = [90, 90];

% Pack into structs.
signals = struct('vw_mps', vw_mps, ...
    'beta_deg', beta_deg, ...
    'phi_deg', phi_deg, ...
    'chi_deg', chi_deg, ...
    'chi_unwrapped_deg', chi_unwrapped_deg);


% Initial values.
init.winch.Lt_m = 1000;
init.winch.w_radps = 0;

end