function Memit = compute_me(C0,Q,V,A,L,D,qn,An,t)
    %     Computes the cumulative mass emitted from the material to air, as
    %     defined in [1] eq 15.
    %
    %   REFERENCES
    %   [1] Deng, B. and C. N. Kim (2004). "An analytical model for VOCs 
    %   emission from dry building materials." Atmospheric Environment 38(8): 1173-1180.
    
    alpha = Q * L^2 / (D * V) ;
    
    qn = qn(:) ;
    An = An(:) ;
    
    comp1 = (((alpha-qn.^2).*sin(qn))./(An.*qn))';
    comp2 = 1 - exp( -D * L^(-2) * qn.^2 * t ) ;
    
    Memit = 2 * C0 * L * A * comp1 * comp2;
    
end

