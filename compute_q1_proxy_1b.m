function q1 = compute_q1_proxy_1b(L,D,K,hm_corrected)
% Calculates q1_proxy using explicit equations, as defined in [1] by eq 17.

% References
% [1] Huang, L. and O. Jolliet (2016). "A parsimonious model for the release
% of chemicals encapsulated in products." Atmospheric Environment 127: 223-235.

    x = hm_corrected*L/(D*K);
    
    Q = 38.7332; B = 1.3386; M = -6.1214; v = 0.0023;
    b = 1/(1+Q*exp(-B*(x-M)))^(1/v);
    a = 1 - b;
            
    y = a*x^0.5 + b*atan(x/(pi/2.2091));

    q1 = y/L;
   
end

