function AIAA_formatting(f, w_fraction, h_fraction)
% Unofficial AIAA formatting made by myself (Jesse).

% AIAA line/textwidth = 469.75502pt
linewidth = 468;  % bit of margin to avoid overfull (somehow).
fontsize(f, 8, 'points')  % Default font size is 10 pt, but should be at least 8 in figures.
grid on

f.Units = 'points';
f.Position = [linewidth, linewidth, w_fraction*linewidth, h_fraction*linewidth];  % [left bottom width height]
end
