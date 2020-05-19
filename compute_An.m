function An = compute_An(Q,V,A,L,D,K,hm, qn )
    %
    %   COMPUTE_AN
    %
    %     Compute An using the complex solution, as defined in [1] by eq 14. 
    %     The ouput is a vector of values for An of same dimension as qn.
    %
    %   REFERENCES
    %   [1] Predicting emissions of volatile and semivolatile organic 
    %       compounds from building materials: A review; Zhe Liu, Wei Yea,
    %       John C. Little; Building and Environment 64 (2013) 7-25.
    %
    
    alpha = Q * L^2 / (D * V) ;
    beta  = A * L / V ;
    Bim    = hm * L / D ;
    
    qn_sq = qn .* qn ;
    
    An = (K*beta + (alpha - qn_sq)*K/Bim + 2) .* qn_sq .* cos( qn ) ...
         + qn .* sin( qn ) .* (K*beta + (alpha - 3*qn_sq)*K/Bim ...
                               + alpha - qn_sq) ;

end
