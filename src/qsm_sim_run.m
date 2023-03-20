%% Load the parameters.
clear
close all
addpath('helper')
addpath('verification_scenarios')
addpath('../parameters/')
[kite, tether, winch, environment] = load_params_mat("my_MegAWES", "..\\parameters");


%% Load a scenario
% TODO: It would be nice if you could also still run the sim from Simulink
% itself without this simIn stuff. Find out how to do that. Do I need a
% similar structure to MegAWES where you first put all variables in your
% workspace and then put then in a simIn object?
% simIn is for temporarily changing the variables of a model. So good to
% use structure of MegAWES.
% simIn = simIn_scenario_1_change_in_wind_speed();

% New method.
% [scenario_name, signals, init, sim_params] = scenario_1A_change_in_path();
[scenario_name, signals, init, sim_params] = scenario_1B_change_in_wind_speed();
% [scenario_name, signals, init, sim_params] = scenario_1C_faster_changes();
% [scenario_name, signals, init, sim_params] = scenario_2_use_massless();


%% Plot path.
close all
figure(1)
plot(signals.vw_mps)
grid on
saveas(gcf, "../results/verification/"+scenario_name+"_vw_mps.png")

plot_path = false;
if plot_path
    figure(2)
    plot(signals.phi_deg.Data, signals.beta_deg.Data);
    grid on
    hold on
    idx = mod(signals.phi_deg.Time, 2) == 0;
    quiver(signals.phi_deg.Data(idx), signals.beta_deg.Data(idx), cosd(signals.angle_deg(idx)), sind(signals.angle_deg(idx)))
    hold off
    axis('equal')
    xlabel('azimuth [deg]')
    ylabel('elevation [deg]')
    legend('path', 'course angle [deg]')
    saveas(gcf, "../results/verification/"+scenario_name+"_path.png")
    
    figure(3)
    plot(signals.chi_deg)
    hold on
    plot(signals.chi_unwrapped_deg)
    grid on
    hold off
end

%% Run sim
% out = sim(simIn);

% New method
tic
out = sim('qsm_sim.slx');
save("../results/verification/" + scenario_name)
toc

%% Analysis
close all

%% tether force over time.
figure(1)
plot(out.Ft_N)
hold on
plot(out.Ft_ref_N, '--')
plot(out.Fte_N)
hold off
legend('actual', 'ref', 'error')
ylabel('tether force [N]')
xlabel('time [s]')
grid on
saveas(gcf, "../results/verification/"+scenario_name+"_tether_force.png")

%% Torque
figure(2)
plot(out.tau_Nm)
hold on
plot(out.tau_PI_Nm)
plot(out.tau_ff_Nm)
legend('combined', 'feedback', 'feedforward')
hold off
ylabel('winch torque [Nm]')
xlabel('time [s]')
grid on
saveas(gcf, "../results/verification/"+scenario_name+"_torque.png")

%% Reeling speed
figure(3)
plot(out.vr_mps)
hold off
ylabel('reeling speed [m/s]')
xlabel('time [s]')
grid on
saveas(gcf, "../results/verification/"+scenario_name+"_reeling_speed.png")

%% Power
figure(4)
plot(out.Pmech_W)
hold on
Pmech_avg_W = mean(out.Pmech_W);
yline(Pmech_avg_W, '--')
hold off
legend('power', 'mean')
ylabel('power output [W]')
xlabel('time [s]')
grid on
saveas(gcf, "../results/verification/"+scenario_name+"_power.png")


%% vr-Ft plot.
figure(5)
plot(out.vr_mps.Data, resample(out.Ft_N, out.vr_mps.Time).Data);
hold on
theta = [146429.5 28806.3 28242.9];
vr_mps = min(0, min(out.vr_mps)):0.1:max(out.vr_mps);
Ft_star_N = theta(1) + theta(2) * vr_mps + theta(3) * vr_mps.^2;
Ft_star_N = max(Ft_star_N, 0.25e6);
Ft_star_N(vr_mps < 0) = 0.25e6;
plot(vr_mps, Ft_star_N, '--');
legend('actual', 'ideal')
hold off
ylabel('tether force [N]')
xlabel('reeling speed [m/s]')
grid on
saveas(gcf, "../results/verification/"+scenario_name+"_vrFt.png")



%% Compare power of two scenario's.
clear
close all
clc


s1 = "scenario_1B_change_in_wind_speed";
s2 = "scenario_2_use_massless";

for s = [s1, s2]
    load("../results/verification/" + s)
    figure(1)
    plot(out.vw_mps.Data, resample(out.Pmech_W, out.vw_mps.Time).Data);
    hold on

    figure(2)
    plot(out.vr_mps.Data, resample(out.Ft_N, out.vr_mps.Time).Data);
    hold on

    figure(3)
    plot(out.Pmech_W)
    Pmech_avg_W = mean(out.Pmech_W)
    yline(Pmech_avg_W, '--')
    hold on
end
figure(1)
legend('my winch', 'massless')
figure(2)
legend('my winch', 'massless')
figure(3)
legend('my winch', 'mean', 'massless', 'mean')


