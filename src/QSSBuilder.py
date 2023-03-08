import os
import sys
from itertools import product

import numpy as np
import pandas as pd
import yaml
from tqdm import tqdm

# When working with submodules it is a bit awkward to import them, especially if
# they're not python packages (with __init__.py files).
# https://stackoverflow.com/questions/29746958/how-to-import-python-file-from-git-submodule
(parent_folder_path, _) = os.path.split(os.getcwd())
workshop_path = os.path.join(parent_folder_path, "workshop")
sys.path.append(workshop_path)
from qsm import Environment, KiteKinematics, SteadyState, SystemProperties
from helper.load_params import load_params


# Class that builds a DataFrame containing lots of possible steady states for
# a certain set of parameters.
# TODO: Path stuff: You now need to run this from the `src` directory and from a Windows machine.
# TODO: logging instead of print, see: https://stackoverflow.com/questions/15727420/using-logging-in-multiple-modules
# and https://docs.python.org/3/howto/logging.html#configuring-logging-for-a-library
class QSSBuilder:
    def __init__(self, param_name: str) -> None:
        self.param_name = param_name
        self.kite, self.tether, _, self.environment = load_params(
            param_name, "..\\parameters"
        )
        self.df = pd.DataFrame()

    @classmethod
    def from_csv(cls, param_name: str):
        # Initialize normally, and then replace the empty DataFrame with the
        # existing one.
        self = cls(param_name)
        self.df = load_df(param_name)
        return self

    def __str__(self) -> str:
        return f"quasi-steady state builder for {self.param_name}"

    def add_steady_states(self, ranges: dict) -> None:
        # Required keys in ranges:
        # vw_mps, Lt_m, phi_deg, beta_deg, chi_deg, Ftk_N for windspeed, tether length,
        # azimuth, elevation, course angle, and kite tether force (control setting).

        # TODO: Make this more flexible: maybe someone wants to put mass in the range.
        #  This could be done by making ranges ALL variables (also wing size etc) and
        #  from this see which one have one value (params) and which ones have multiple
        #  (ranges) and then assign appropriately. Then also need two savefiles one with
        #  fixed params and one with the DataFrame. Can't all put it in the DataFrame
        #  because it will get too big.

        # `range` is a dictionary with iterables where the keys will be the column names
        # and the rows will be all possible combinations of the ranges of those keys
        df_new = self.make_df_new(ranges)

        # Merge the new DataFrame into `self.df`, where if SteadyStates were already
        # calculated for a certain state it is preserved and doesn't need to be
        # calculated again.
        self.df_new_into_self_df(df_new)

        # For the rows where 'calculates_SS' is false we must now calculate the
        # SteadyState.
        self.calculate_new_states()

    @staticmethod
    def make_df_new(ranges: dict, verbose=True) -> pd.DataFrame:
        if verbose:
            print(
                "Making new DataFrame with all possible combinations in `ranges`.",
                end="\r",
            )

        df_new = pd.DataFrame(product(*ranges.values()), columns=ranges.keys())

        if verbose:
            print("Done making new DataFrame with all possible combinations in ranges.")
        return df_new

    def df_new_into_self_df(self, df_new: pd.DataFrame, verbose=True) -> None:
        if verbose:
            print(f"Merging new values into existing DataFrame", end="\r")
            total_rows_before_merge = self.df.shape[0]

        # Combine the new and existing DataFrame.
        if self.df.empty:
            # If the existing one is empty we just put the `df_new` in its place.
            self.df = df_new
            self.df["calculated_SS"] = False
        else:
            # Merge the new DataFrame with all new ranges into the existing one, merge
            # allows duplicates to be skipped. When a state in `df_new` is new, it will
            # get NaN in all other columns (calculated_SS, f, Ftg_N, P_W). Then replace
            # the NaN value in calculated_SS with `False` so that we know that we have
            # to run the steady state calculation on it.
            # TODO: Maybe concatenating the dataframes and then removing dupicates is
            # faster?
            self.df = self.df.merge(df_new, how="outer", on=list(df_new.columns))
            self.df["calculated_SS"].fillna(False, inplace=True)

        if verbose:
            print(
                f"Merged {self.df.shape[0] - total_rows_before_merge} new rows into \
                  the DataFrame."
            )

    def calculate_new_states(self) -> None:
        # tqdm shows a nice progress bar while we loop over the entire DataFrame.
        total_rows = self.df.shape[0]
        for row in tqdm(
            self.df.itertuples(),
            desc="Calculating quasi-steady states",
            total=total_rows,
        ):

            # Don't recalculate previously calculated steady states.
            if row.calculated_SS:
                continue

            # Calculate reeling factor, tether force ground, and power.
            f, Ftg_N, P_W, E_eff = self.calculate_new_state(row)

            # Put values in the DataFrame.
            self.df.at[row.Index, "calculated_SS"] = True
            self.df.at[row.Index, "f"] = f
            self.df.at[row.Index, "Ftg_N"] = Ftg_N
            self.df.at[row.Index, "P_W"] = P_W
            self.df.at[row.Index, "E_eff"] = E_eff
            self.df.at[row.Index, "vr_mps"] = f * row.vw_mps

    def calculate_new_state(self, row):
        # Wrapper around SteadyState class of qsm.py
        # Make the environment, system properties and kite kinematic objects.
        env_state = {
            "wind_speed": row.vw_mps,
            "air_density": self.environment["rho_kgpm3"],
        }
        env_state = Environment(**env_state)

        sys_props = {
            "kite_projected_area": self.kite["S_m2"],
            "kite_mass": self.kite["m_kg"],
            "tether_density": self.tether["rho_kgpm3"],
            "tether_diameter": self.tether["r_m"] * 2,
            "kite_lift_coefficient_powered": self.kite["C_L"],
            "kite_drag_coefficient_powered": self.kite["C_L"] / self.kite["E"],
            "tether_drag_coefficient": self.tether["Cd_t"],
        }
        sys_props = SystemProperties(sys_props)

        kite_kinematics = {
            "straight_tether_length": row.Lt_m,
            "azimuth_angle": np.deg2rad(row.phi_deg),
            "elevation_angle": np.deg2rad(row.beta_deg),
            "course_angle": np.deg2rad(row.chi_deg),
        }
        kite_kinematics = KiteKinematics(**kite_kinematics)

        # The mass of the tether and aerodynamic properties including the tether length
        # can now be calculated.
        sys_props.update(kite_kinematics.straight_tether_length)

        ss = SteadyState()
        # I found that the control setting tether force kite was most stable.
        ss.control_settings = ("tether_force_kite", row.Ftk_N)

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

        # It sometimes happens that the SteadyState converges but returns None.
        if ss.tether_force_kite is None:
            result = (np.nan, np.nan, np.nan, np.nan)
            ss.converged = False

        if ss.converged:
            # Check whether the state it found was actually the state we wanted.
            # Especially when the control setting was the ground tether force this
            # didn't always work due to some issues in qsm.py.
            # Even np.isclose sometime failed because something was not right with the
            # input types.
            try:
                np.isclose(
                    ss.tether_force_kite, row.Ftk_N, rtol=1e-5
                ), f"Error in `qsm.py`: Steady state converged but final tether force \
                 ({ss.tether_force_kite}) is not equal to the setpoint ({row.Ftk_N}).)"
            except Exception as e:
                print(e)
                print(
                    f"final tether force ({ss.tether_force_kite}, setpoint ({row.Ftk_N})."
                )
                result = (np.nan, np.nan, np.nan, np.nan)

        return result

    def save_df(self, name: str = None, verbose=True) -> None:
        if not name:
            name = self.param_name

        if verbose:
            print(f"Saving DataFrame of {self} as {name}.csv")

        self.df.to_csv(f"../results/{name}.csv")


def load_df(param_name, verbose=True):
    if verbose:
        print(f"Loading {param_name}.csv", end="\r")

    if os.path.isfile(f"../results/{param_name}.csv"):
        df = pd.read_csv(f"../results/{param_name}.csv", index_col=0)
    else:
        raise FileNotFoundError(f"Can't find {param_name}.csv in results folder.")

    if verbose:
        print(f"Done loading {param_name}.csv")

    return df
