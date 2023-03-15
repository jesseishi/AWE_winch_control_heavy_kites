function E_eff = calc_E_eff(tether_length, kite, tether)
% calc_E_Eff: calculate the effective lift-to-drag ratio of a kite, taking
% the tether drag into account.
CD = kite.CL / kite.E;
CD_eff = CD + 0.5 * tether.r_m * tether_length / kite.S_m2 * tether.Cd_t;
E_eff = kite.CL / CD_eff;
end
