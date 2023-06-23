function [scenario_name, signals, init, sim_params] = scenario_2A_power_curve()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
scenario_name = mfilename;

% One loop = 20 sec.
t_end = 1000;
sim_params = struct('t_end', t_end, ...
    'use_massless_controller', false, ...
    'use_larger_winch_PI', false, ...
    'use_zgraggen', false, ...
    'use_massless_kitemodel', false);

% Define a very approximate figure of eight movement.
t = (0:0.1:t_end)';
% Beta starts going down so the kite flies up on the outer part of the
% figure-of-eight.
beta_deg = 30 - 10 * sin(pi/10*t);
phi_deg = 0 + 40 * sin(pi/10/2*t);  % azimuth is halve the frequency
% We can calculate the angle by taking atan2(Delta_phi, Delta_beta).
% When then needs a transformation to chi which starts downward and goes
% clockwise.
angle_rad = atan2(-pi*cos(pi/10*t), 2*pi * cos(pi/10/2*t));

% angle_deg is discontinuous at -180 and 180. For interpolation to work
% properly when this signal is imported to Simulink, it needs to be
% continuous.
angle_rad = unwrap(angle_rad);
angle_deg = rad2deg(angle_rad);
chi_unwrapped_deg = 270 - angle_deg;
chi_deg = wrapTo360(chi_unwrapped_deg);

% In timeseries
beta_deg = timeseries(beta_deg, t);
phi_deg = timeseries(phi_deg, t);
chi_unwrapped_deg = timeseries(chi_unwrapped_deg, t);
chi_deg = timeseries(chi_deg, t);

% Need quite a high windspeed because the quasi-steady model cannot handle
% accelerations which is what would happen at low wind speed while flying
% down.
vw_mps = linspace(10, 30, length(t))';
vw_mps = timeseries(vw_mps, t);

% Pack into structs.
signals = struct('vw_mps', vw_mps, ...
    'beta_deg', beta_deg, ...
    'phi_deg', phi_deg, ...
    'chi_deg', chi_deg, ...
    'chi_unwrapped_deg', chi_unwrapped_deg, ...
    'angle_deg', angle_deg);


% Initial states.
init.winch.Lt_m = 1000;
init.winch.w_radps = 3.5 / 1.5;

end