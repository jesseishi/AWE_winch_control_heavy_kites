# AWE_winch_control_heavy_kites
Analytic analysis of how to make optimal winch controller for kites where mass cannot be neglected. Used for my MSc thesis. The goal is to gain insight into how a kite with mass required different winch control than a kite without mass, for which the theory is readily available. Optimal control techniques can also be used for this goal but in this work a more analytic approach is taken, by using the quasi-steady model from [1].

## References:
- [1] R. van der Vlugt, A. Bley, M. Noom, and R. Schmehl, “Quasi-steady model of a pumping kite power system,” Renewable Energy, vol. 131, pp. 83–99, Feb. 2019, doi: 10.1016/j.renene.2018.07.023.

## Structure
- parameters
  - Folder with different `.yaml` files that define the parameters of different systems that can be analysed.
- results
  - Folder to store results.
- src
  - `qsm_sim.slx`: Simulink model using the quasi-steady model (qsm) model from `workshop`. It calls the python file `qsm_wrapper_for_simulink.py` from Simulink to use the quasi-steady model.
  - `qsm_wrapper_for_simulink.py`: see `qsm_sim.slx`.
  - `qss_heavy_analysis.ipynb`: jupyter notebook that investigates the quasi-steady states for a AWE system with mass. It has a focus on creating optimal winch control strategies. It uses the `QSSBuilder` class to generate all the data for a lot of different steady-states.
  - `qss_massless_analysis.ipynb`: Same as `qss_heavy_analysis.ipynb` but for kites without mass.
  - `QSSBuilder_test.ipynb`: test file for `QSSBuilder.py`.
  - `QSSBuilder.py`: Defines the `QSSBuilder` class which can loop through a lot of different steady-state conditions and call the quasi-steady state model from `workshop` to calculate the steady-state. It then stores the results in a Pandas DataFrame and in the results folder as a `.csv` file.
  - `winch_parameter_analysis.m`: Analysis the effect of different winch parameters.
  - `winch_PI_force_controller.m`: Design of a winch PI force controller.
- workshop: Code from https://github.com/awecourse/workshop, imported as submodule. This code implements the quasi-steady model, as described by https://doi.org/10.1016/j.renene.2018.07.023.

## Tips
- I always had my IDE open at the base folder. It might be that some relative paths need to be adjusted if you place it elsewhere. Except for SimpleSim.slx which wouldn't work that way (which has something to do with calling Python but I have no idea what - just glad that it works now).