%% Set the parameters.
clear
close all
addpath('helper')
[kite, tether, winch, environment] = load_params_mat("my_MegAWES", "../parameters");


%% Update some parameters, based on tether length.
Lt_m = 1500;  % tether length -> for worst-case.
kite.E_eff = calc_E_eff(Lt_m, kite, tether);
kite.CR_eff = kite.CL * sqrt(1 + 1/kite.E_eff^2);
kite.C = 0.5 * environment.rho_kgpm3 * kite.S_m2 * kite.CR_eff * (1 + kite.E_eff^2);


%% Make a state-space model with inputs torque and wind speed.
% There are two state-space models, one with a massless model and one with a linear approximate model of a heavy kite.
use_massless = true;
if use_massless

    % Trim condition.
    Ft_0 = 0.5e6;
    beta_0 = deg2rad(30);
    phi_0 = deg2rad(35);
    f_0 = 1/3 * cos(beta_0) * cos(phi_0);
    vw_0 = sqrt(Ft_0 / (kite.C * (cos(beta_0)*cos(phi_0) - f_0)^2));
    vr_0 = vw_0 * f_0;

    % Assert that that calculation is corrent.
    assert(abs(Ft_0 - kite.C * (vw_0*cos(beta_0)*cos(phi_0) - vr_0)^2) < 0.001)
    assert(abs(Ft_0 - kite.C * vw_0^2*(cos(beta_0)*cos(phi_0) - f_0)^2) < 0.001)

    % Temporary variable as this expression comes back often.
    temp = -2 * kite.C * (vw_0*cos(beta_0)*cos(phi_0) - vr_0);

    % Build the state space system.
    A = [(-winch.friction + winch.r_m^2*temp)/winch.J_kgm2, 0;
         -temp,                                             0];
    B_tau = [-winch.r_m / winch.J_kgm2;
             0];
    B_vw = [-winch.r_m^2/winch.J_kgm2*temp*cos(beta_0)*cos(phi_0);
            temp*cos(beta_0)*cos(phi_0)];
    B = [B_tau, B_vw];
    C = [A(2, :);
         0, 1   ];
    D_tau = [B_tau(2, :);
             0 ];
    D_vw = [B_vw(2, :);
            0];
    D = [D_tau, D_vw];
    sys = ss(A, B, C, D, 'StateName', {'v_r', 'Fte_int'}, 'InputName', {'\tau', 'v_w'}, 'OutputName', {'F_t_e', 'Fte_int'});
else
    theta_1 = -290436.9;  % massless: -113 100
    theta_2 =  210452.0;  % massless:   80 236
    A = [(theta_1*winch.r_m^2 - winch.friction)/winch.J_kgm2, 0;
         -theta_1,                                            0];
    B_tau = [-winch.r_m/winch.J_kgm2;
             0];
    B_vw = [theta_2*winch.r_m^2/winch.J_kgm2;
            -theta_2];
    B = [B_tau, B_vw];
    C = [-theta_1, 0;
         0,        1];
    D_tau = [0;
             0];
    D_vw = [-theta_2;
            0];
    D = [D_tau, D_vw];
    sys = ss(A, B, C, D, 'StateName', {'v_r', 'Fte_int'}, 'InputName', {'\tau', 'v_w'}, 'OutputName', {'F_t_e', 'Fte_int'});
end


%% Analysis
close all
figure(1)
step(sys(1, 2))  % From wind speed to tether force. It is already stable, but we want it faster and without steady-state error.
grid('on')
saveas(gcf, "../results/open_loop_step_massless_"+use_massless+".png")

figure(2)
bode(sys(1, 2))
grid('on')
xlim([1e-1, 1e3])
saveas(gcf, "../results/open_loop_bode_massless_"+use_massless+".png")

%% PID design.
Ft_over_tau = sys(1, 1);
pidTuner(Ft_over_tau)
% Here we want output (Ft) disturbance rejection. I don't know how to
% select the bandwidth appropriately.
% Selected: bandwidth: 30 rad/s, phase margin 80 deg.

%% PID parameters
if use_massless
    Kp = -1.481;
    Ki = -53.53;
else
    Kp = -0.41768;
    Ki = -47.9043;
end
% For tau = - K y:
K = [Kp, Ki];

%% Verify
% The closed loop is now Ft output with input Ft error.
s = tf('s');
Ft_over_Fte = feedback(Ft_over_tau*(Kp + Ki/s), 1);
Ft_over_Fte.InputName = 'Ftref';
figure(4)
step(Ft_over_Fte)
hold on
step(feedback(Ft_over_tau, 1));
grid on
legend('closed loop', 'open loop')
hold off
saveas(gcf, "../results/closed_loop_step_ftref_massless_"+use_massless+".png")

%% Closed loop.
% Closed loop response with disturbance vw (input tau doesn't exist anymore
% as it is now used for feedback).
sys_cl_vw = ss(A - B_tau*K*C, B_vw - B_tau*K*D_vw, C, D_vw, 'StateName', {'v_r', 'Fte_int'}, 'InputName', 'v_w', 'OutputName', {'F_t_e', 'Ft_int'});

%% Analysis
% Check the response when there is a step in windspeed.
figure(5)
step(sys_cl_vw(1))
hold on
step(sys(1, 2))
legend('closed loop', 'open loop')
hold off
grid('on')
saveas(gcf, "../results/closed_loop_step_massless_"+use_massless+".png")

% And to a ramp in wind.
t = linspace(0, 100, 10000);
% u = t';
% u = t.^2;
% u = sin(0.01*2*pi*t.^2);  % Poor man's bode plot.
u = sin(pi/10*t);  % w_MegAWES = pi/10 rad/s
% u = u + 0.01*randn(size(u));  % turbulence, but since there is not
% damping in the tether and there's a feedforward term D this does not
% produce realistic results.
% TODO: open loop should be also good here right?
figure(6)
lsim(sys_cl_vw(1), u, t)
hold on
lsim(sys(1, 2), u, t);
legend('closed loop', 'open loop')
hold off
grid('on')
saveas(gcf, "../results/closed_loop_at_MegAWES_freq_massless_"+use_massless+".png")

% Bode plots
figure(7)
bode(sys_cl_vw(1))
hold on
bode(sys(1, 2))
legend('closed loop', 'open loop')
hold off
grid('on')
xlim([1e-1, 1e3])
saveas(gcf, "../results/closed_loop_bode_massless_"+use_massless+".png")

% Step in tether force reference.




