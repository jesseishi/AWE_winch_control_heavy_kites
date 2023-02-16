import numpy as np

from QSSBuilder import QSSBuilder

qss = QSSBuilder.from_csv("MegAWES_SixDofPaper")

qss.add_steady_states(
    {
        "vw_mps": np.unique(np.linspace(0, 30, 13, dtype=int)),
        "Lt_m": np.unique(np.linspace(500, 1500, 3, dtype=int)),
        "phi_deg": np.unique(np.linspace(0, 30, 3, dtype=int)),
        "beta_deg": np.unique(np.linspace(0, 45, 4, dtype=int)),
        "chi_deg": np.unique(np.linspace(0, 180, 5, dtype=int)),
        # TODO: Smarter way to find different Ftk_N that work than looping over all values.
        # Use dtype int for the merge into an existing DataFrame to be less error-prone.
        # Otherwise it sometimes couldn't match two floats with each other. Integer
        # precision is also more than enough here.
        # "Ftk_N": np.unique(np.logspace(1, 7, 103, dtype=int)),
        "Ftk_N": np.unique(np.arange(10, 1e7, 50000, dtype=int)),
    }
)

qss.save_df()
