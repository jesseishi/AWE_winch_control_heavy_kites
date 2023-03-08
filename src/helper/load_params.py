import os
import yaml


def load_params(param_name, params_dir=None):

    if params_dir is None:
        path = f"{param_name}.yaml"
    else:
        path = os.path.join(params_dir, f"{param_name}.yaml")

    # Load the parameter file. These are values that will remain fixed.
    with open(path, "r") as f:
        try:
            params = yaml.safe_load(f)
            kite = params["kite"]
            tether = params["tether"]
            winch = params["winch"]
            environment = params["environment"]

        except yaml.YAMLError as exc:
            print(exc)

    return kite, tether, winch, environment


# For testing.
if __name__ == "__main__":
    param_name = "MegAWES"
    kite, tether, winch, environment = load_params(param_name, "parameters")
    print(kite, tether, winch, environment)
