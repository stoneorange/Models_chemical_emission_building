function y_t = compute_y( C0, Q,V,A,L,D, hm, qn, An, t )
    %
    %   COMPUTE_Y
    %
    %     Compute approximation (partial sum) of y(t) as defined in [1] 
    %   by eq 13. The number of terms in the partial sum is given by the
    %   number of elements of qn and An.
    %
    %   NOTE
    %     The sum is obtained by matrix multiplication between two
    %   components of its argument:
    %       comp1 = qn .* sin(qn) ./ An
    %       comp2 = exp( -D * L^(-2) * qn.^2 * t )
    %   with comp1 a 1 x n_qn vector, and qn.^2 * t and hence comp2 a 
    %   n_qn x n_t matrix.    
    %
    %   REFERENCES
    %   [1] Predicting emissions of volatile and semivolatile organic 
    %       compounds from building materials: A review; Zhe Liu, Wei Yea,
    %       John C. Little; Building and Environment 64 (2013) 7-25.
    %
    
    alpha = Q * L^2 / (D * V) ;
    beta  = A * L / V ;
    Bim    = hm * L / D ;
    
    % - Enforce qn and An to be column vectors, and t to be a row vector.
    qn = qn(:) ;
    An = An(:) ;
    t  = t(:).' ;
    
    % - Compute comp1 and comp2.
    comp1 = (qn .* sin( qn ) ./ An).' ;
    comp2 = exp( -D * L^(-2) * qn.^2 * t ) ;

    % - Compute y_t.
    y_t = 2 * C0 * beta * comp1 * comp2 ;
    
end