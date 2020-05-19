function C_xt = compute_c( C0, Q,V,A,L,D, hm, qn, An, t, x )
    %     Computes the concentration in the material at depth x and time t
    %     using the complex solution, as defined in [1] by eq 12.
    %
    %   REFERENCES
    %   [1] Predicting emissions of volatile and semivolatile organic 
    %       compounds from building materials: A review; Zhe Liu, Wei Yea,
    %       John C. Little; Building and Environment 64 (2013) 7-25.
    
    alpha = Q * L^2 / (D * V) ;
    beta  = A * L / V ;
    Bim    = hm * L / D ;


    qn = qn(:) ;
    An = An(:) ;
    n  = length(t);
    
    comp3 = (alpha-qn.^2)./An;
    comp3_2 = repmat (comp3, 1, n);
    comp4 = exp( -D * L^(-2) * qn.^2 * t ) ;
    comp5 = cos((x/L)* qn');
    
    C_xt = 2 * C0 * comp5 * (comp3_2.*comp4);
    

end

