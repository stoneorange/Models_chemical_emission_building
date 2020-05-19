function Memit = compute_me(C0,Q,V,A,L,D,qn,An,t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    alpha = Q * L^2 / (D * V) ;
    
    qn = qn(:) ;
    An = An(:) ;
    
    comp1 = (((alpha-qn.^2).*sin(qn))./(An.*qn))';
    comp2 = 1 - exp( -D * L^(-2) * qn.^2 * t ) ;
    
    Memit = 2 * C0 * L * A * comp1 * comp2;
    
end

