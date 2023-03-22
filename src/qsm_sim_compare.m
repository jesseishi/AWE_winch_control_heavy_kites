% Compare the power output of my winch controller and PI tuning vs what you
% would get with a winch controller and PI tuning based on massless flight.

clear
close all

s1 = "scenario_2A_power_curve";
s2 = "scenario_2B_power_curve_massless_controller";

for s = [s1, s2]
    load("../results/verification/" + s)
    figure(1)
    idx = mod(out.vw_mps.Time, 1) == 0;
    plot(out.vw_mps.Data, resample(out.Pmech_W, out.vw_mps.Time).Data);
    t = 0:0.1:1000;
    hold on
    plot(resample(out.vw_mps, t).Data, movmean(resample(out.Pmech_W, t).Data, 200), '--')

    figure(2)
    plot(out.vr_mps.Data, resample(out.Ft_N, out.vr_mps.Time).Data);
    hold on

    figure(3)
    plot(out.Pmech_W)
    Pmech_avg_W = mean(out.Pmech_W, 'Weighting', 'Time')  % The datapoints are not evenly spread out over time, so weight by time to correct.
%     yline(Pmech_avg_W, '--')
    hold on
    plot(t, movmean(resample(out.Pmech_W, t).Data, 200), '--')
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
legend('heavy strategy', 'moving mean', 'massless', 'moving mean')
grid on
xlabel('time [s]')
ylabel('power output [W]')
saveas(gcf, "../results/verification/compare_to_massless_power_output.png")
grid on




