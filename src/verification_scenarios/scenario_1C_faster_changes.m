function [scenario_name, signals, init, sim_params] = scenario_1C_faster_changes()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
scenario_name = mfilename;

% One loop = 20 sec.
t_end = 5;
sim_params = struct('t_end', t_end, ...
    'use_massless_controller', false, ...
    'use_larger_winch_PI', false, ...
    'use_zgraggen', false, ...
    'use_massless_kitemodel', false);

% Set the windspeed.
% Fastest wind oscillation is 10*pi rad/s -> period = 0.2
t = (0:0.002:t_end)';
vw_mps = 15 + 15*sin(pi/10*t) + 5*sin(pi*t) + 4*sin(5*pi*t);
vw_mps = 15 + 1*sin(10*pi*t);
% But never below 0.
vw_mps = max(vw_mps, 0.001);  % QSM model needs at least a little bit of wind.
vw_mps = timeseries(vw_mps, t, 'Name', 'wind speed [m/s]');

% Other variables.
% These variables are variable in other verification scenario's so they
% must be compatible with `from_workspace`, meaning we must do this awkward
% construction of [constant, constant].
beta_deg = [0, 0];
phi_deg = [0, 0];
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
