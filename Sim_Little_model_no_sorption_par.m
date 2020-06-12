function [ Mf_emit_end, C_air] = Sim_Little_model_no_sorption_par( hm,Am,dm,Km,Q,Kpa,TSP,C0,t )
% Simplify the single diffusion model by assuming chemical concentration in
% building material is constant, and no sorption on other indoor surfaces.
% Adapted from Little et al., 2012.

% Ref: Little, J. C., et al. (2012). "Rapid methods to estimate potential 
% exposure to semivolatile organic compounds in the indoor environment." 
% Environmental Science & Technology 46(20): 11171-11178.

C_air = Am*hm*C0/Km/(Am*hm/(1+Kpa*TSP)+Q);

Mf_emit_end = Am*hm*(C0/Km - C_air/(1+Kpa*TSP))*t / (Am*dm*C0);


end
