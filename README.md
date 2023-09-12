# AWE_winch_control_heavy_kites
This repository was used for my MSc thesis [1]. In this repository I worked towards two goals:
1. Finding a method for sizing the winch of a ground-gen airborne wind energy system. More specifically, to find an upper bound on the winch sizing constant (inertia / radius^2).
2. Analytic analysis of how to make optimal winch control curve for kites where mass cannot be neglected. This was done by analysing the optimal reel-out speed and tether force for kites with mass, by using the quasi-steady model from [2].

## Structure
- parameters
  - Folder with different `.yaml` files that define the parameters of different systems that can be analysed.
- results
  - Folder to store results.
- src
  - `helper` folder: bunch of helper functions for other files.
  - `verification_scenarios`: folder with different scenario's to run on the `qsm_sim`.
  - `presentation_plots.ipynb`: jupyter notebook to make an interactive plot to find the optimal winch control curve for kites with mass.
  - `qsm_sim_compare.m`: compare two or more completed runs from the `qsm_sim`.
  - `qsm_sim_run.m`: run the `qsm_sim` with a certain scenario and analyse the results.
  - `qsm_sim_slx`: Simulink model using the quasi-steady model (qsm) model from `workshop` which implements the model from [2].
  - `qsm_wrapper_for_simulink.py`: since the quasi-steady model with mass from [2] is written in python, this file provides a wrapper function around the model that can be called from Simulink.
  - `qss_heavy_analysis.ipynb`: jupyter notebook that investigates the quasi-steady states for a AWE system with mass. It has a focus on creating optimal winch control strategies. It uses the `QSSBuilder` class to generate all the data for a lot of different steady-states.
  - `qss_massless_analysis.ipynb`: Same as `qss_heavy_analysis.ipynb` but for kites without mass.
  - `QSSBuilder_test.ipynb`: test file for `QSSBuilder.py`.
  - `QSSBuilder.py`: Defines the `QSSBuilder` class which can loop through a lot of different steady-state conditions and call the quasi-steady state model from `workshop` to calculate the steady-state. It then stores the results in a Pandas DataFrame and in the results folder as a `.csv` file.
  - `winch_parameter_analysis.m`: Analysis of the effect of the size of the winch on an airborne wind energy system. This can be used to derive an upper limit on the winch sizing constant (inertia / radius^2).
  - `winch_PI_force_controller.m`: Design of a winch PI force controller.
- workshop: Code from https://github.com/awecourse/workshop, imported as submodule. This code implements the quasi-steady model with mass from [2].

## References:
- [1] J. Hummel, “Kite Tether Force Control: Reducing Power Fluctuations for Utility-Scale Airborne Wind Energy Systems,” Delft University of Technology, Delft, Netherlands, 2023. [Online]. Available: http://resolver.tudelft.nl/uuid:2a0b5cb1-ab4d-44bf-b6b6-929c612941fb
- [2] R. van der Vlugt, A. Bley, M. Noom, and R. Schmehl, “Quasi-steady model of a pumping kite power system,” Renewable Energy, vol. 131, pp. 83–99, Feb. 2019, doi: 10.1016/j.renene.2018.07.023.
