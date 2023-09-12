%% Set some winch parameters.
clear
close all
addpath('helper')
[kite, tether, winch, environment] = load_params_mat("MegAWES", "../parameters");
s = tf('s');

%% Update some parameters, based on tether length.
Lt_m = 1000;  % tether length.
kite.E_eff = calc_E_eff(Lt_m, kite, tether);
kite.CR_eff = kite.CL * sqrt(1 + 1/kite.E_eff^2);
kite.C = 0.5 * environment.rho_kgpm3 * kite.S_m2 * kite.CR_eff * (1 + kite.E_eff^2);


%% Trim point analysis.
vr = 0:0.001:5;
Ft = (4*kite.C*winch.r_m^2*vr.^2+winch.friction*vr ) / winch.r_m^2;
figure
plot(vr, Ft)


%% Bode plot of P(s)/F(s) for different values of J.
% Uses Ft as input (Ft is thus not affected by vr -> would need a kite
% model for that and now all effects on the kite (wind, beta, phi, chi,
% Lt, etc..) are collapsed into the effect of Ft. The winch does assume a
% massless kite model to provide feedforward control.
figure
hold on

vr_0 = 4.0278932;  % 0.5 MN tether force.
% vr_0 = 2.8487939;  % 0.25 MN tether force.
Ft_0 = (4*kite.C*winch.r_m^2*vr_0^2+winch.friction*vr_0 ) / winch.r_m^2

for J = logspace(3, 7, 5)
    tf_PF = (12*kite.C*vr_0^2*winch.r_m^2) / (J*s + 8*kite.C*winch.r_m^2*vr_0 + winch.friction);

    % P and F are not directly comparable so it's better to make a plot of
    % P(s)/P_ideal / F(s).
    tf_PF_normalized = (12*kite.C*vr_0^2*winch.r_m^2) / (J*s + 8*kite.C*winch.r_m^2*vr_0 + winch.friction) * ...
        (8*kite.C*winch.r_m^2*vr_0 + winch.friction) / (12*kite.C*vr_0^2*winch.r_m^2);
    bode(tf_PF_normalized)
%     damp(tf_PF_normalized);
end
grid on
xline(pi/10)
legend('J = 1e3 kgm^2', 'J = 1e4 kgm^2', 'J = 1e5 kgm^2', 'J = 1e6 kgm^2', 'J = 1e7 kgm^2', 'MegAWES frequency')

% TODO: idk how to put the xline in both upper and lower plots
% automatically (now I just click and redo the command..).
% saveas(gcf, '../results/bode_P_over_Pideal_over_f.png')


%% Visualise the magnitude of the bode plot at the MegAWES frequency for different J, r and v_r_0.
w = pi/10;
% w = w*10;  % Safety factor so that we can also respond to e.g. gusts. TODO: Do we want this? Maybe it's better to not respond to that. Maybe a safety factor of 2 would be better.


N = 67;
Jv = logspace(1, 7, N);
% v_r_0v = linspace(0.1, 20, N);
v_r_0v = logspace(log10(0.01), log10(20), N);  % logspace for better resolution at low numbers.
rv = linspace(0.1, 5, N);
[JV, RV, V_R_0V] = meshgrid(Jv, rv, v_r_0v);

MAGS_V = abs( (8*kite.C*RV.^2.*V_R_0V + winch.friction) ./ (JV*w*1i + 8*kite.C*RV.^2.*V_R_0V + winch.friction));

figure('Units', 'centimeters', 'Position', [5, 5, 15.75, 10])
v_r_slice = [0.1, 1, 2.848, 4.0278932, 20];
colors = {'red', 'blue', 'cyan', 'green', 'magenta'};
for i=1:length(v_r_slice)
    % The order of [0.95, 0.99] doesn't matter, it will sort it from low to
    % high. This is unfortunate because we want to display it the other way
    % around, but whatever.
    p = contourslice(JV, RV, V_R_0V, MAGS_V, [], [], v_r_slice(i), [sqrt(2)/2, 0.95, 0.99]);  % sqrt(2)/2 is the magnitude of the bandwidth frequency.
    p(1).LineStyle = ':';
    p(2).LineStyle = '--';
    p(3).LineStyle = '-';
    p(1).EdgeColor = colors{i};
    p(2).EdgeColor = colors{i};
    p(3).EdgeColor = colors{i};
    p(1).DisplayName = "v_r_0 = " + v_r_slice(i) + ", pf = " + p(1).UserData;
    p(2).DisplayName = "v_r_0 = " + v_r_slice(i) + ", pf = " + p(2).UserData;
    p(3).DisplayName = "v_r_0 = " + v_r_slice(i) + ", pf = " + p(3).UserData;
end

legend('Location', 'eastoutside')
xlabel('inertia [kgm^2]'); ylabel('radius [m]')
grid on
set(gca, 'xscale', 'log')

xlim([Jv(1), Jv(end)])

% title('inertia and radius for certain power fractions (pf) at distinct reel-out speeds.')
saveas(gcf, '../results/power_frac_Jrvr0.png')

%% Chosen values:
J = 1e4;
r = 1.5;
vr_0 = 1;

pole = -(8 * kite.C * r^2 * vr_0 + winch.friction) / J;

P_norm_over_Ft = (12*kite.C*vr_0^2*r^2) / (J*s + 8*kite.C*r^2*vr_0 + winch.friction) * ...
        (8*kite.C*r^2*vr_0 + winch.friction) / (12*kite.C*vr_0^2*r^2);
damp(P_norm_over_Ft)

figure(1)
bodemag(P_norm_over_Ft)
grid on
xline(pi/10, '--', 'MegAWES frequency')
xline(-pole, '--', 'bandwidth')
yline(mag2db(0.99), '--', '0.99 power fraction')
ylim([-20, 1])
saveas(gcf, '../results/bode_P_over_Pideal_over_f_tuned_winch.png')

% % Create zoomed plot inside the figure.
% axes('Position',[.55 .452 .2 .2])
% box on
figure(2)
bodemag(P_norm_over_Ft, logspace(log10(pi/10/2), log10(pi/10*2), 10))
xlim([pi/10/2, pi/10*2])
ylim([mag2db(0.99)*2, 0])
% title('zoom')
xline(pi/10, '--', 'MegAWES frequency')
xline(-pole, '--', 'bandwidth')
yline(mag2db(0.99), '--', '0.99 power fraction')
% xlabel('')
% ylabel('')
grid on
saveas(gcf, '../results/bode_P_over_Pideal_over_f_tuned_winch_zoom.png')
