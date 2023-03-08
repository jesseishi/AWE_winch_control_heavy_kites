# Used by Matlab to access the qsm.py library.
import os
import sys

import numpy as np

from load_params import load_params

# When working with submodules it is a bit awkward to import them, especially if
# they're not python packages (with __init__.py files).
# https://stackoverflow.com/questions/29746958/how-to-import-python-file-from-git-submodule
(parent_folder_path, _) = os.path.split(os.getcwd())
workshop_path = os.path.join(parent_folder_path, "workshop")
sys.path.append(workshop_path)
from qsm import Environment, KiteKinematics, SteadyState, SystemProperties


def calc_Ft(
    vw_mps,
    vr_mps,
    Lt_m,
    beta_deg,
    phi_deg,
    chi_deg,
    params_dir=None,
):
    kite, tether, _, environment = load_params("my_MegAWES", params_dir)

    env_state = {
        "wind_speed": vw_mps,
        "air_density": environment["rho_kgpm3"],
    }
    env_state = Environment(**env_state)

    sys_props = {
        "kite_projected_area": kite["S_m2"],
        "kite_mass": kite["m_kg"],
        "tether_density": tether["rho_kgpm3"],
        "tether_diameter": tether["r_m"] * 2,
        "kite_lift_coefficient_powered": kite["C_L"],
        "kite_drag_coefficient_powered": kite["C_L"] / kite["E"],
        "tether_drag_coefficient": tether["Cd_t"],
    }
    sys_props = SystemProperties(sys_props)

    kite_kinematics = {
        "straight_tether_length": Lt_m,
        "azimuth_angle": np.deg2rad(phi_deg),
        "elevation_angle": np.deg2rad(beta_deg),
        "course_angle": np.deg2rad(chi_deg),
    }
    kite_kinematics = KiteKinematics(**kite_kinematics)

    # The mass of the tether and aerodynamic properties including the tether length
    # can now be calculated.
    sys_props.update(kite_kinematics.straight_tether_length)

    ss = SteadyState()
    # I found that the control setting tether force kite was most stable.
    ss.control_settings = ("reeling_speed", vr_mps)

    try:
        ss.find_state(sys_props, env_state, kite_kinematics)
        result = (
            ss.reeling_factor,
            ss.tether_force_ground,
            ss.power_ground,
            sys_props.lift_to_drag,
        )
    except Exception:
        result = (np.nan, np.nan, np.nan, np.nan)

    return result[1]


if __name__ == "__main__":
    print(calc_Ft(10, 1, 1000, 0, 0, 0, params_dir="parameters"))
