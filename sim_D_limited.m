function [ Mf_emit, C_gas, C_surf, C_gas_avg ] = sim_D_limited(Dm, Km, dm, Am, C0, hm, Q, t)
% The D-limited model, as defined in [1].

% References
% [1] Huang, L. and O. Jolliet (2016). "A parsimonious model for the release
% of chemicals encapsulated in products." Atmospheric Environment 127: 223-235.

    hm_corrected = hm/(1+Am*hm/Q);  % Eq 13 in [1]
    q1_proxy = compute_q1_proxy_1b(dm,Dm,Km,hm_corrected);  % Eq 17 in [1]
    M0 = C0 * Am * dm; 

    % Table 1 in [1]
    X = q1_proxy^2 * Dm;
    alpha = max(1 - 1483*X^0.45023, 0.85);
    beta1 = max(0.01366*X^(-0.20891)*q1_proxy, q1_proxy);
    beta2 = exp(0.0091*log(X)^2 + 0.2115*log(X) + 2.591) * beta1;
    
    % Eq 20 in [1]
    Mf_emit = alpha*(1 - exp(-beta1^2*Dm*t)) + (1-alpha)*(1 - exp(-beta2^2*Dm*t));
    
    % C_gas is the concentration in gas phase of indoor air.
    % C_gas = dMf_emit/dt * M0 / Q_cor.
    C_gas = (alpha*beta1^2*Dm*exp(-beta1^2*Dm*t) + (1-alpha)*beta2^2*Dm*exp(-beta2^2*Dm*t)) * M0 / Q;
    
    % C_surf is the concentration at the material surface.
    % C_surf = dMf_emit/dt * M0 * Km/(Am*hm_corrected) Eq 21 in [1]
    C_surf = (alpha*beta1^2*Dm*exp(-beta1^2*Dm*t) + (1-alpha)*beta2^2*Dm*exp(-beta2^2*Dm*t)) * C0*dm*Km / hm_corrected;
    
    % C_gas_avg = integrate(C_gas(t)dt)/V/T = Mf_emit * M0 / Q / T.
    C_gas_avg = Mf_emit * C0*Am*dm / t / Q;
    
end

