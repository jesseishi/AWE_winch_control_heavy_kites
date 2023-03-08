function E_eff = calc_E_eff(tether_length, kite, tether)
% calc_E_Eff: calculate the effective lift-to-drag ratio of a kite, taking
% the tether drag into account.
CD = kite.C_L / kite.E;
CD_eff = CD + 0.5 * tether.r_m * tether_length / kite.S_m2 * tether.Cd_t;
E_eff = kite.C_L / CD_eff;
end
