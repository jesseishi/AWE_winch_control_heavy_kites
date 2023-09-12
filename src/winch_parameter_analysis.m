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


%% Trim points analysis.
vr = 0:0.001:9;
Ft = 4*kite.C*vr.^2;
figure
plot(vr, Ft)
xlabel('Reel-out speed, m/s')
ylabel('Tether force, N')
ylim([0, Ft(end)])
xlim([0, 9])
grid on

if 1
    hold on

    Ft_max = 1.39e6;
    P_max = Ft_max * sqrt(Ft_max / (4*kite.C));
    Ft_P_max = P_max ./ vr;
    
%     yline(Ft_max, '-', 'F_t max')
    plot(vr, Ft_max*ones(size(vr)))
    plot(vr, Ft_P_max)
    legend('winch control curve', 'tether force limit', 'power limit')

    ylim([0, Ft(end)])
end

AIAA_formatting(gcf, 0.6, 0.6/1.6)

% saveas(gcf, '../results/paper/vrFt.png')
saveas(gcf, '../results/paper/vrFt_with_limits.png')


%% Pole position analysis.
beta_0 = 0;
phi_0 = 0;
vr_0 = 5;
vr_over_vw = tf(4*kite.C*winch.r_m^2*vr_0*cos(beta_0)*cos(phi_0), ...
    [winch.J_kgm2, 12*kite.C*winch.r_m^2*vr_0 + winch.friction]);
damp(vr_over_vw)
sprintf('pole is at %.1e', -(12*kite.C*winch.r_m^2*vr_0 + winch.friction)/winch.J_kgm2)

%% Bode plot of P(s)/F(s) for different values of J.
% Uses Ft as input (Ft is thus not affected by vr -> would need a kite
% model for that and now all effects on the kite (wind, beta, phi, chi,
% Lt, etc..) are collapsed into the effect of Ft. The winch does assume a
% massless kite model to provide feedforward control.
figure
hold on

vr_0 = 5;

for J = logspace(3, 7, 5)
    tf_PF = (48*kite.C^2*vr_0^3*winch.r_m^2*cos(beta_0)*cos(phi_0)) / ...
        (J*s + 12*kite.C*winch.r_m^2*vr_0 + winch.friction);

    % P and F are not directly comparable so it's better to make a plot of
    % P(s)/P_ideal / F(s).
    P_ideal = (48*kite.C^2*vr_0^3*winch.r_m^2*cos(beta_0)*cos(phi_0)) / ...
        (12*kite.C*winch.r_m^2*vr_0 + winch.friction);
    tf_PF_normalized = tf_PF / P_ideal;
    bode(tf_PF_normalized)
%     damp(tf_PF_normalized);
end
grid on
xline(pi/10, '--')
legend('J = 1e3 kgm^2', 'J = 1e4 kgm^2', 'J = 1e5 kgm^2', 'J = 1e6 kgm^2', 'J = 1e7 kgm^2', 'MegAWES frequency')

% TODO: idk how to put the xline in both upper and lower plots
% automatically (now I just click and redo the command..).
% saveas(gcf, '../results/bode_P_over_Pideal_over_vw.png')


%% Visualise the magnitude of the bode plot at the MegAWES frequency for different J, r and v_r_0.
w = pi/10;
w = w*10;  % Safety factor so that we can also respond to e.g. gusts. TODO: Do we want this? Maybe it's better to not respond to that. Maybe a safety factor of 2 would be better.


N = 67;  % Ensures all 10^N where N is an integer are in the logspace for Jv.
Jv = logspace(1, 7, N);
v_r_0v = logspace(log10(0.01), log10(20), N);  % logspace for better resolution at low numbers.
rv = linspace(0.5, 5, N);
[JV, RV, V_R_0V] = meshgrid(Jv, rv, v_r_0v);

MAGS_PF = abs ((48*kite.C^2*V_R_0V.^3.*RV.^2*cos(beta_0)*cos(phi_0)) ./ ...
        (JV*1i*w + 12*kite.C.*RV.^2.*V_R_0V + winch.friction));

MAGS_P_ideal = (48*kite.C^2*V_R_0V.^3.*RV.^2*cos(beta_0)*cos(phi_0)) ./ ...
        (12*kite.C.*RV.^2.*V_R_0V + winch.friction);

MAGS_V = MAGS_PF ./ MAGS_P_ideal;

% figure('Units', 'centimeters', 'Position', [5, 5, 15.75, 10])
figure('Units', 'centimeters', 'Position', [5, 5, 15.75, 15.75])
% v_r_slice = [0.1, 1, 2.848, 4.0278932, 20];
v_r_slice = [1, 4, 10];
colors = {'red', 'blue', 'cyan', 'green', 'magenta'};
for i=1:length(v_r_slice)
    % The order of [0.95, 0.99] doesn't matter, it will sort it from low to
    % high. This is unfortunate because we want to display it the other way
    % around, but whatever.
    p = contourslice(JV, RV, V_R_0V, MAGS_V, [], [], v_r_slice(i), [sqrt(2)/2, 0.95, 0.99]);  % sqrt(2)/2 is the magnitude of the bandwidth frequency.
%     p = contourslice(JV, RV, V_R_0V, MAGS_V, [], [], v_r_slice(i), [0.95, 0.99]);  % sqrt(2)/2 is the magnitude of the bandwidth frequency.
    p(1).LineStyle = ':';
    p(2).LineStyle = '--';
    p(3).LineStyle = '-';
    p(1).EdgeColor = colors{i};
    p(2).EdgeColor = colors{i};
    p(3).EdgeColor = colors{i};
    p(1).DisplayName = sprintf('reel-out speed: %4.1f m/s, power fraction: %.2f', v_r_slice(i), p(1).UserData);
    p(2).DisplayName = sprintf('reel-out speed: %4.1f m/s, power fraction: %.2f', v_r_slice(i), p(2).UserData);
    p(3).DisplayName = sprintf('reel-out speed: %4.1f m/s, power fraction: %.2f', v_r_slice(i), p(3).UserData);
end

legend('Location', 'southoutside')
xlabel('inertia [kgm^2]'); ylabel('radius [m]')
grid on
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')

xlim([Jv(1), Jv(end)])
ylim([rv(1), rv(end)])
yticks([1:5])

% title('inertia and radius for certain power fractions (pf) at distinct reel-out speeds.')
saveas(gcf, '../results/power_frac_Jrvr0.png')

%% What if the inertia was 5.2e7 kgm2 (such as that one 5MW wind turbine)
w = pi/10;
w = w*10;  % Safety factor so that we can also respond to e.g. gusts. TODO: Do we want this? Maybe it's better to not respond to that. Maybe a safety factor of 2 would be better.


N = 67;
Jv = logspace(1, 8, N);
v_r_0v = logspace(log10(0.01), log10(20), N);  % logspace for better resolution at low numbers.
rv = linspace(0.1, 150, N);
[JV, RV, V_R_0V] = meshgrid(Jv, rv, v_r_0v);

MAGS_PF = abs ((48*kite.C^2*V_R_0V.^3.*RV.^2*cos(beta_0)*cos(phi_0)) ./ ...
        (JV*1i*w + 12*kite.C.*RV.^2.*V_R_0V + winch.friction));

MAGS_P_ideal = (48*kite.C^2*V_R_0V.^3.*RV.^2*cos(beta_0)*cos(phi_0)) ./ ...
        (12*kite.C.*RV.^2.*V_R_0V + winch.friction);

MAGS_V = MAGS_PF ./ MAGS_P_ideal;

% figure('Units', 'centimeters', 'Position', [5, 5, 15.75, 10])
figure('Units', 'centimeters', 'Position', [5, 5, 15.75, 15.75])
% v_r_slice = [0.1, 1, 2.848, 4.0278932, 20];
v_r_slice = [0.1, 1, 5, 20];
colors = {'red', 'blue', 'cyan', 'green', 'magenta'};
for i=1:length(v_r_slice)
    % The order of [0.95, 0.99] doesn't matter, it will sort it from low to
    % high. This is unfortunate because we want to display it the other way
    % around, but whatever.
%     p = contourslice(JV, RV, V_R_0V, MAGS_V, [], [], v_r_slice(i), [sqrt(2)/2, 0.95, 0.99]);  % sqrt(2)/2 is the magnitude of the bandwidth frequency.
    p = contourslice(JV, RV, V_R_0V, MAGS_V, [], [], v_r_slice(i), [0.95, 0.99]);  % sqrt(2)/2 is the magnitude of the bandwidth frequency.
    p(1).LineStyle = '--';
    p(2).LineStyle = '-';
    p(1).EdgeColor = colors{i};
    p(2).EdgeColor = colors{i};
    p(1).DisplayName = sprintf('reel-out speed: %4.1f m/s, power fraction: %.2f', v_r_slice(i), p(1).UserData);
    p(2).DisplayName = sprintf('reel-out speed: %4.1f m/s, power fraction: %.2f', v_r_slice(i), p(1).UserData);
end

legend('Location', 'southoutside')
xlabel('inertia [kgm^2]'); ylabel('radius [m]')
grid on
set(gca, 'xscale', 'log')

xlim([Jv(1), Jv(end)])


%% Chosen values:
vr_0 = 1.0;
J_kgm2 = 1e4;
r_m = 1.5;

tf_PF = (48*kite.C^2*vr_0^3*r_m^2*cos(beta_0)*cos(phi_0)) / ...
    (J_kgm2*s + 12*kite.C*r_m^2*vr_0 + winch.friction);
P_ideal = (48*kite.C^2*vr_0^3*r_m^2*cos(beta_0)*cos(phi_0)) / ...
    (12*kite.C*r_m^2*vr_0 + winch.friction);
tf_PF_normalized = tf_PF / P_ideal;
damp(tf_PF_normalized)
bandwidth = (12*kite.C*r_m^2*vr_0 + winch.friction)/J_kgm2;

figure(1)
bodemag(tf_PF_normalized)
grid on
xline(pi/10, '--', 'MegAWES frequency', 'LabelVerticalAlignment', 'bottom')
xline(pi/10*10, '--', 'MegAWES frequency with safety margin', 'LabelVerticalAlignment', 'bottom')
xline(bandwidth, '--', 'bandwidth', 'LabelVerticalAlignment', 'bottom')
yline(mag2db(0.99), '--', '0.99 power fraction')
yline(0, '-k')
ylim([-5, 0.5])
xlim([pi/15, 2*bandwidth])
saveas(gcf, '../results/bode_P_over_Pideal_over_vw_tuned_winch.png')


%% For paper: better illustrations.
% I can introduce a winch sizing constant: Kw = J/r^2 and ignore friction.
w = pi/10;
w = w*10;

N = 1e3;
vrs = logspace(log10(0.01), log10(10), N);
Kws = logspace(3, 6, N);

[VRS, KWS] = meshgrid(vrs, Kws);

MAGS = abs(12 * kite.C * VRS ./ (KWS * 1j* w + 12 * kite.C * VRS));

figure()
[M, c] = contourf(VRS, KWS, MAGS, linspace(0, 1, 101));  % a bit slow..
c.EdgeColor = 'none';
hold on
contour(VRS, KWS, MAGS, [0.90, 0.95, 0.99], '-k', 'ShowText', true, "LabelFormat", "%.2f")

colormap(parula)
clim([0.5, 1.0])  % This can really mess with the plot. It doesn't
% necessarily adjust the colors of the contourf.
c = colorbar;
c.Label.String = 'Normalized power output, -';
xlabel('Trim reel-out speed, m/s')
ylabel('Winch sizing constant, kg')
% set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
AIAA_formatting(gcf, 0.5, 0.5/1.6)

saveas(gcf, '../results/paper/winch_sizing.png')


%% Bode plot for vr/vw
r_m = 1.5;
J_kgm2 = 1e4;

vr0 = 1;

vr_over_vw = tf(4 * kite.C * r_m^2 * vr0, ...
    [J_kgm2, 12 * kite.C * r_m^2 * vr0 + winch.friction]);

P_norm = vr_over_vw * 3;

figure(1)
bodemag(P_norm)
ylabel('Magnitude, db')
xlabel('Frequency, rad/s')
AIAA_formatting(gcf, 0.6, 0.6/1.2)
saveas(gcf, '../results/paper/winch_bode.png')

%% Winch control strategy.
vr = linspace(0, 7, 100);
Ft_star = 4 * kite.C * vr.^2;

figure;
plot(vr, Ft_star)
xlabel('Reel-out speed, m/s')
ylabel('Tether force, N')
AIAA_formatting(gcf, 0.6, 0.6/1.6)
saveas(gcf, '../results/paper/vrFt.png')

%% Tether force as function of vw.
% THIS REALLY ADDS A LOT TO THIS ANALYSIS.
r_m = 0.4;
J_kgm2 = 32;

vr0 = 5;
vw0 = 3*vr0;  % f = 1/3.

vr_over_vw = tf(4 * kite.C * r_m^2 * vr0, ...
    [J_kgm2, 12 * kite.C * r_m^2 * vr0 + winch.friction]);

D = 2*kite.C*(vw0 - vr0);
Ft_over_vw = D * (1-vr_over_vw);

% Use simplified normalization - valid for small friction.
tf_Ft = Ft_over_vw/(D*2/3);
tf_Pn = 3*vr_over_vw;

close all
figure(1)
bode(tf_Ft, tf_Pn)
legend('F_t', 'P_n')
AIAA_formatting(gcf, 0.6, 0.6)
saveas(gcf, '../results/paper/winch_bode2.png')

figure(2)
step(tf_Ft, tf_Pn)
legend('F_t', 'P_n')
AIAA_formatting(gcf, 0.6, 0.6)
saveas(gcf, '../results/paper/winch_step2.png')


%% MAGS map but for Ft.
% I can introduce a winch sizing constant: Kw = J/r^2 and ignore friction.
w = pi/10;
w = w*10;

N = 1e3;
vrs = logspace(log10(0.01), log10(10), N);
Kws = logspace(3, 6, N);

[VRS, KWS] = meshgrid(vrs, Kws);

% MAGS = abs(12 * kite.C * VRS ./ (KWS * 1j* w + 12 * kite.C * VRS));
% MAGS = abs((4*kite.C*VRS.*KWS*1j*w + 32*kite.C^2*r_m^2*VRS.^2) ./ ...
%     (KWS*1j*w + 12*kite.C*VRS));  % Magnitude of tether force TF.
% MAGS = MAGS ./ ((32*kite.C^2*r_m^2*VRS.^2) ./ 12*kite.C*VRS);  % Normalise.

MAGS = abs(1.5 * (KWS*1j*w + 8*kite.C*VRS) ./ (KWS*1j*w + 12*kite.C*VRS));

% figure('Units', 'centimeters', 'Position', [5, 5, 15.75, 15.75])
figure()
% Nice to have both continuous surface and then lines on the contours but
% since surface is 3D you don't have a grid...
% s = surface(VRS, KWS, MAGS);
% set(s, 'EdgeColor', 'none')
% hold on
% contour3(VRS, KWS, MAGS, [db2mag(-3), 0.90, 0.95, 0.99], '-k');

% contourf(VRS, KWS, MAGS, [0.5, db2mag(-3), 0.90, 0.95, 0.99]);

[M, c] = contourf(VRS, KWS, MAGS, linspace(1, 1.5, 101));  % a bit slow..
c.EdgeColor = 'none';
hold on
contour(VRS, KWS, MAGS, [1.01, 1.05, 1.10], '-k', 'ShowText', true, "LabelFormat", "%.2f")

colormap(flipud(parula))
clim([1.0, 1.5])  % This can really mess with the plot. It doesn't
% necessarily adjust the colors of the contourf.
c = colorbar;
c.Label.String = 'Normalized tether force, -';
xlabel('Trim reel-out speed, m/s')
ylabel('Winch sizing constant, kg')
% set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
AIAA_formatting(gcf, 0.5, 0.5/1.6)

saveas(gcf, '../results/paper/winch_sizing_Ft.png')


%% Nice visual for in a presentation. Can also output gifs to explain the effect of lag caused by the inertia.
figure
LineWidth = 3;

vr = 0:0.1:6;
Ft = 4*kite.C*vr.^2;
plot(vr, Ft, '--k', 'LineWidth', LineWidth);
hold on
grid on
xlabel('Reel-out speed, m/s')
ylabel('Tether Force, N')
ylim([-2e5, 12e5])
legend('ideal', 'Location', 'NorthWest')

saveas(gcf, '../results/presentation/winch_sizing_0.png')

r_m = 1.5;

t = 0:0.01:30;
base_frequency_rad = pi/10;

rJf = [1.5, 1e6, 1;
       0.1, 1e4, 1;
%        0.1, 1e4, 10;
       1.5, 1e4, 1];
legendnames = {'too large inertia', '', '', 'too low radius', '', '', 'final parameters', '', ''};

colors = get(gca,'ColorOrder');

for i = 1:size(rJf, 1)
    r_m = rJf(i, 1);
    J_kgm2 = rJf(i, 2);
    freq_multiplier = rJf(i, 3);

    for vr0 = 1:2:5
    
        % Trim condition.
        vw0 = 3*vr0;  % f = 1/3.
        Ft0 = 4*kite.C*vr0^2;
        
        vr_over_vw = tf(4 * kite.C * r_m^2 * vr0, ...
            [J_kgm2, 12 * kite.C * r_m^2 * vr0 + winch.friction]);
        D = 2*kite.C*(vw0 - vr0);
        Ft_over_vw = D * (1-vr_over_vw);
    
        
%         u = 1.5 * (sin(freq_multiplier*0.05*2*pi*t - pi/2) + 1);  % Phase offset and bias make the sine start at 0 and oscilate between 0 and 2.
        u = 2*sin(freq_multiplier*0.05*2*pi*t);

        vr = lsim(vr_over_vw, u, t);
        Ft = lsim(Ft_over_vw, u, t);

        % Only when making gifs:
%         legend off
%         for ti = 1:length(t)
%             if mod(ti, 100) == 0
%                 plot(vr(1:ti) + vr0, Ft(1:ti) + Ft0, 'Color', colors(i, :), 'LineWidth', LineWidth)
%                 exportgraphics(gcf, sprintf('windDisturbance_%d.gif', i), 'Append', true);
%             end
%         end
    
        plot(vr + vr0, Ft + Ft0, 'Color', colors(i, :), 'LineWidth', LineWidth)
    end
    legend(['ideal', legendnames(1:i*3)])
    saveas(gcf, sprintf('../results/presentation/winch_sizing2_%d.png', i))
end


% AIAA_formatting(gcf, 0.8, 0.8/1.6)

%% Another plot for the presentation.
figure
hold on
grid on
xlabel('Reel-out speed, m/s')
ylabel('Tether Force, N')
% ylim([-2e5, 12e5])
% legend('ideal', 'Location', 'NorthWest')

vr = 0:0.01:9.0;
Ft = 4*kite.C*vr.^2;
P = vr.*Ft;

mymap = [0.4660 0.6740 0.1880
         0.9290 0.6940 0.1250
         0.8500 0.3250 0.0980];
P_idx = 1e6*[0, 10, 15];
for i = 1:3
    idx = P>=P_idx(i);
    scatter(vr(idx), Ft(idx), [], mymap(i, :), 'filled')
end

legend('  0 MW < P < 10 MW', '10 MW < P < 15 MW', '15 MW < P < 20 MW', ...
    'Location', 'NorthWest')

saveas(gcf, '../results/presentation/power_lim.png')

Ft_P10MW = (10e6 * sqrt(4 * kite.C))^(2/3);
yline(Ft_P10MW)

legend('  0 MW < P < 10 MW', '10 MW < P < 15 MW', '15 MW < P < 20 MW', ...
    'Tether Force Limit', ...
    'Location', 'NorthWest')

saveas(gcf, '../results/presentation/power_force_lim.png')

%% Overview plot.
figure
hold on
grid on
xlabel('Reel-out speed, m/s')
ylabel('Tether Force, N')

plot(vr, Ft, 'LineWidth', 3)
legend('without force limit')
saveas(gcf, '../results/presentation/vrFt_without_limit.png')
Ft_limited = min(Ft, Ft_P10MW);
plot(vr, Ft_limited, 'LineWidth', 3)
legend('without force limit', 'with force limit')
saveas(gcf, '../results/presentation/vrFt_with_limit.png')


%% Overview plot 2.
vr = linspace(0, 9, 100);
Ft = linspace(0, 2.5e6, 100);
[VR, FT] = meshgrid(vr, Ft);

P = VR .* FT;
figure
Ft = 4*kite.C*vr.^2;
Ft_limited = min(Ft, Ft_P10MW);
plot(vr, Ft_limited)  % Cheating color coding.
hold on
plot(vr, Ft_limited, 'LineWidth',3)
grid on

% surface(VR, FT, P, 'EdgeColor','none')
contour(VR, FT, P./1e6, [0, 5, 10, 15, 20],"ShowText",true,"LabelFormat","%0.1f MW",'LineWidth',3)
legend('','Optimal with tether force limit', '')
xlabel('Reel-out speed, m/s')
ylabel('Tether force, N')

saveas(gcf, '../results/presentation/vrFt_with_powerlimit.png')
