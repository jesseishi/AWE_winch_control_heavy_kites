function [kite, tether, winch, environment] = load_params_mat(param_name, param_dir)
% Load a yaml file with kite, tether, winch and environment parameters.
% To keep the parameters easy to read and adjust, they are stored in yaml
% files. However, Matlab can't really read those so I'm calling python from
% here to read them for us.
pyrun("from helper.load_params import load_params");
[kite, tether, winch, environment] = pyrun(...
    "kite, tether, winch, environment = load_params(param_name, param_dir)", ...
    ["kite", "tether", "winch", "environment"], ...
    param_name=param_name, param_dir=param_dir);

% Turn Python dicts into structs.
kite = struct(kite);
tether = struct(tether);
winch = struct(winch);
environment = struct(environment);
end
