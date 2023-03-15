%% Set the parameters.
clear
close all
addpath('src')
addpath('src/helper')
[kite, tether, winch, environment] = load_params_mat("my_MegAWES", "parameters");

%% Update some paramters, based on tether length.
Lt_m = 1000;  % tether length.
kite.E_eff = calc_E_eff(Lt_m, kite, tether);
kite.CR_eff = kite.C_L * sqrt(1 + 1/kite.E_eff^2);
kite.C = 0.5 * environment.rho_kgpm3 * kite.S_m2 * kite.CR_eff * (1 + kite.E_eff^2);



%% Make a state-space model with inputs torque and wind speed.
% There are two state-space models, one with a massless model and one with a linear approximate model of a heavy kite.
use_massless = false;
if use_massless
    % Trim conditions.
    Ft_0 = 0.25e6;
    vw_0 = sqrt(9*Ft_0/(4*kite.C));
    vr_0 = vw_0/3;
    
    % State-space system with extra state Ft_int which can be used for
    % feedback.
    A = [((-2*vw_0 + 2*vr_0) * kite.C * winch.r^2 - winch.friction)/winch.J, 0;
         kite.C * (-2*vw_0 + 2*vr_0),                                        0];
    % Split inputs of torque and wind because one is an input and the otherone
    % is a disturbance.
    B_tau = [-winch.r / winch.J;
             0                 ];
    B_vw = [((2*vw_0 - 2*vr_0) * kite.C * winch.r^2) / winch.J;
            kite.C * (2*vw_0 - 2*vr_0)];
    B = [B_tau, B_vw];
    C = [kite.C * (-2*vw_0 + 2*vr_0), 0;
         0,                             1];
    D_tau = [0;
             0];
    D_vw = [kite.C * (2*vw_0 - 2*vr_0);
            0];
    D = [D_tau, D_vw];
    
    sys = ss(A, B, C, D, 'StateName', {'vr', 'Ft_int'}, 'InputName', {'tau', 'vw'}, 'OutputName', {'Ft', 'Ft_int'});
    sys_tau = ss(A, B_tau, C, D_tau, 'StateName', {'vr', 'Ft_int'}, 'InputName', 'tau', 'OutputName', {'Ft', 'Ft_int'});
else
    theta_1 = -150842.6;
    theta_2 = 109489.2;
    A = [(theta_1*winch.r^2 - winch.friction)/winch.J, 0;
         -theta_1,                                     0];
    B_tau = [-winch.r/winch.J;
             0];
    B_vw = [theta_2*winch.r^2/winch.J;
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
step(sys(1, 2))  % From wind speed to tether force. It is already stable, but we want it faster.
grid('on')
saveas(gcf, 'results/open_loop_step.png')

figure(2)
bode(sys(1, 2))
grid('on')
saveas(gcf, 'results/open_loop_bode.png')

%% PID design.
Ft_over_tau = sys(1, 1);
pidTuner(Ft_over_tau)
% Here we want output (Ft) disturbance rejection at a frequency slighly
% lower than 31.415 (frequency of system with feedforward control law).

%% PID parameters
Kp = -0.9;
Ki = -30;
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
legend('closed loop', 'open loop')
hold off

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
saveas(gcf, 'results/closed_loop_step.png')

% And to a ramp in wind.
t = linspace(0, 100, 10000);
% u = t';
% u = t.^2;
% u = sin(0.01*2*pi*t.^2);  % Poor man's bode plot.
u = sin(pi/10*t);  % w_MegAWES = pi/10 rad/s
% u = u + 0.01*randn(size(u));  % turbulence, but since there is not
% damping in the tether and there's a feedforward term D this does not
% produce realistic results.
figure(6)
lsim(sys_cl_vw(1), u, t)
hold on
lsim(sys(1, 2), u, t);
legend('closed loop', 'open loop')
hold off
grid('on')
saveas(gcf, 'results/closed_loop_at_MegAWES_freq.png')

% Bode plots
figure(7)
bode(sys_cl_vw(1))
hold on
bode(sys(1, 2))
legend('closed loop', 'open loop')
hold off
grid('on')
saveas(gcf, 'results/closed_loop_bode.png')



