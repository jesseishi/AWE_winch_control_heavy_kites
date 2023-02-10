# AWE_winch_control_heavy_kites
Analytic analysis of how to make optimal winch controller for kites where mass cannot be neglected. Used for my MSc thesis. The goal is to gain insight into how a kite with mass required different winch control than a kite without mass, for which the theory is readily available. Optimal control techniques can also be used for this goal but in this work a more analytic approach is taken, by using the quasi-steady model from [1].

## References:
- [1] R. van der Vlugt, A. Bley, M. Noom, and R. Schmehl, “Quasi-steady model of a pumping kite power system,” Renewable Energy, vol. 131, pp. 83–99, Feb. 2019, doi: 10.1016/j.renene.2018.07.023.

## Structure
- results
- src
- workshop: Code from https://github.com/awecourse/workshop, imported as submodule. This code implements the quasi-steady model, as described by https://doi.org/10.1016/j.renene.2018.07.023.
