%% Load the parameters.
clear
close all
addpath('helper')
addpath('verification_scenarios')
addpath('../parameters/')
[kite, tether, winch, environment] = load_params_mat("my_MegAWES", "..\\parameters");


%% Load a scenario
[scenario_name, signals, init, sim_params] = scenario_1A_change_in_path();
% [scenario_name, signals, init, sim_params] = scenario_1B_change_in_wind_speed();
% [scenario_name, signals, init, sim_params] = scenario_1C_faster_changes();
% [scenario_name, signals, init, sim_params] = scenario_2A_power_curve();
% [scenario_name, signals, init, sim_params] = scenario_2B_power_curve_massless_controller();
% [scenario_name, signals, init, sim_params] = scenario_2C_constant_wind();
% [scenario_name, signals, init, sim_params] = scenario_2D_constant_wind_massless_controller();
% [scenario_name, signals, init, sim_params, winch] = scenario_3A_bigger_inertia(winch);
% [scenario_name, signals, init, sim_params] = scenario_4A_massless_kitemodel();
% [scenario_name, signals, init, sim_params, winch] = scenario_4B_massless_kitemodel_bigger_inertia(winch);
% [scenario_name, signals, init, sim_params, winch] = scenario_4C_massless_kitemodel_smaller_radius(winch);
% [scenario_name, signals, init, sim_params, winch] = scenario_4D_massless_kitemodel_BIGGER_inertia(winch);


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
% ylim([-0.5e6, 5e6])

if scenario_name == "scenario_1A_change_in_path"
    ylim([-0.5e6, 1.5e6])
    % Create zoomed plot inside the figure.
    axes('Position',[.55 .452 .2 .2])
    box on
    idx = (out.Ft_N.Time > 25.2) & (out.Ft_N.Time < 26.3);
    plot(out.Ft_N.Time(idx), out.Ft_N.Data(idx))
    hold on
    plot(out.Ft_ref_N.Time(idx), out.Ft_ref_N.Data(idx), '--')
    hold off
    ylim([0.49e6, 0.52e6])
    title('zoom')
    grid
elseif scenario_name == "scenario_1B_change_in_wind_speed"
    axes('Position',[.675 .30 .2 .2])
    box on
    idx = (out.Ft_N.Time > 30) & (out.Ft_N.Time < 38);
    plot(out.Ft_N.Time(idx), out.Ft_N.Data(idx))
    hold on
    plot(out.Ft_ref_N.Time(idx), out.Ft_ref_N.Data(idx), '--')
    hold off
    ylim([0.45e6, 0.55e6])
    title('zoom')
    grid
end


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
% ylim([-1e7, 1e7])
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
plot(out.vr_mps * out.Ft_N)
idx = out.Pmech_W.Time > 0.1;  % Calculate the average power after it has reached the steady-state.
temp = timeseries(out.Pmech_W.Data(idx), out.Pmech_W.Time(idx));
Pmech_avg_W = mean(temp, 'Weighting', 'Time');
yline(Pmech_avg_W, '--', sprintf('%.2f MW', Pmech_avg_W/1e6))


hold off
legend('power out', 'power in', 'mean')
ylabel('power output [W]')
xlabel('time [s]')
grid on
ylim([-2e7, 7e7])
saveas(gcf, "../results/verification/"+scenario_name+"_power.png")


%% vr-Ft plot.
figure(5)
plot(out.vr_mps.Data, resample(out.Ft_N, out.vr_mps.Time).Data);
hold on
theta = [187486.8, 7784.0, 30865.2];
vr_mps = min(0, min(out.vr_mps)):0.1:max(out.vr_mps);
Ft_star_N = theta(1) + theta(2) * vr_mps + theta(3) * vr_mps.^2;
Ft_star_N = max(Ft_star_N, 0.5e6);
Ft_star_N(vr_mps < 0) = 0.5e6;
plot(vr_mps, Ft_star_N, '--');

% Only for massless.
% Lt_m = 1000;  % tether length.
% kite.E_eff = calc_E_eff(Lt_m, kite, tether);
% kite.CR_eff = kite.CL * sqrt(1 + 1/kite.E_eff^2);
% kite.C = 0.5 * environment.rho_kgpm3 * kite.S_m2 * kite.CR_eff * (1 + kite.E_eff^2);
% plot(vr_mps, 4*kite.C*vr_mps.^2, '--');
legend('actual', 'ideal')
hold off
ylabel('tether force [N]')
xlabel('reeling speed [m/s]')
grid on
saveas(gcf, "../results/verification/"+scenario_name+"_vrFt.png")


%% Is the system lagging?
figure(6)
plot(out.vr_mps)
hold on
plot(signals.vw_mps)
ylabel('speed [m/s]')

% yyaxis right
% ylabel('reel-out factor [-]')
% plot(out.vr_mps.Time, out.vr_mps.Data ./ resample(signals.vw_mps, out.vr_mps.Time).Data)
% ylim([0, 1])
% yline(1/3, '--', 'optimal reel-out factor')

grid on
xlabel('time [s]')
legend('reel-out speed', 'wind speed', 'reel-out factor')
saveas(gcf, "../results/verification/"+scenario_name+"_vrvw.png")


