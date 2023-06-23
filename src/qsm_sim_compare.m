% Compare the power output of my winch controller and PI tuning vs what you
% would get with a winch controller and PI tuning based on massless flight.

clear
close all

% s1 = "scenario_2A_power_curve";
% s2 = "scenario_2B_power_curve_massless_controller";
s1 = "scenario_2C_constant_wind";
s2 = "scenario_2D_constant_wind_massless_controller";
ss = [s1, s2];

aligns = {'top', 'bottom'};

for i = 1:length(ss)
    s = ss(i);
    align = aligns{i};

    load("../results/verification/" + s)
    figure(1)
    idx = mod(out.vw_mps.Time, 1) == 0;
    plot(out.vw_mps.Data, resample(out.Pmech_W, out.vw_mps.Time).Data);
    t = 0:0.1:40;
    hold on
    plot(resample(out.vw_mps, t).Data, movmean(resample(out.Pmech_W, t).Data, 200), '--')

    figure(2)
    plot(out.vr_mps.Data, resample(out.Ft_N, out.vr_mps.Time).Data);
    hold on

    figure(3)
    plot(out.Pmech_W)
    Pmech_avg_W = mean(out.Pmech_W, 'Weighting', 'Time')  % The datapoints are not evenly spread out over time, so weight by time to correct.
%     yline(Pmech_avg_W, '--', sprintf('%.2f', Pmech_avg_W))
    hold on
%     plot(t, movmean(resample(out.Pmech_W, t).Data, 200), '--')
%     hold on

    figure(4)
    plot(out.vr_mps)
    vr_mean_mps = mean(out.vr_mps, 'Weighting', 'Time');
    yline(vr_mean_mps, '--', sprintf('%.2f', vr_mean_mps), 'LabelVerticalAlignment', align);
    hold on

    figure(5)
    plot(out.vr_mps ./ resample(signals.vw_mps, out.vr_mps.Time))
    hold on
end
figure(1)
legend('my winch', 'my winch', 'massless', 'massless')

figure(2)
legend('heavy strategy', 'massless')
grid on
xlabel('reel-out speed [m/s]')
ylabel('tether force [N]')
saveas(gcf, "../results/verification/compare_to_massless_vrFt.png")

figure(3)
% legend('heavy strategy', 'moving mean', 'massless', 'moving mean')
legend('heavy strategy', 'massless')
grid on
xlabel('time [s]')
ylabel('power output [W]')
saveas(gcf, "../results/verification/compare_to_massless_power_output.png")
grid on

figure(4)
legend('heavy strategy', '', 'massless', '')
grid on
xlabel('time [s]')
ylabel('reeling speed [m/s]')
saveas(gcf, "../results/verification/compare_to_massless_reeling_speed.png")

figure(5)
ylabel('reel-out factor [-]')
grid on
saveas(gcf, "../results/verification/compare_to_massless_reeling_factor.png")
