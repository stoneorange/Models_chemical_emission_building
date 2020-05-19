function [qn, pole, hasCorrectedErrors, internals] = ...
    compute_qn_vector( Q, V, A, L, D, K, hm, nRoots, correctErrors, ...
    bisCntMax, tol )
    %
    %   COMPUTE_QN_VECTOR
    %
    %      Compute qn roots as defined in [1], eq. 18, using an approach
    %   based on a vector bisection.
    %
    %   NOTES
    %     As developed in [2], the function
    %
    %                   f(qn) = qn .* tan(qn)
    %
    %   is continuous and strictly monotonically increasing in every 
    %   interval [-pi/2,pi/2]+k*pi while the function
    %
    %       g(qn) = (alpha-qn.^2) ./ (beta*K + (alpha-qn.^2)*K/Bim) 
    %
    %   is continuous and strictly monotonically decreasing everywhere [3]
    %   except at its unique pole qn_pole. Because f and g have opposite 
    %   monotonicity in each interval [-pi/2,pi/2]+k*pi not containing
    %   qn_pole, and because f(qn) ranges from -inf to +inf there, the
    %   equation
    %                        f(qn) = g(qn)
    %
    %   will have exactly one solution in each of those intervals. These
    %   solutions are identically the roots of f(qn)-g(qn). The anomalous
    %   interval [-pi/2,pi/2]+k0*pi containing the pole qn_pole can be
    %   handled by splitting it into 2 subintervals to the left and right
    %   of qn_pole. Each of those subintervals must have exactly one root
    %   by the same monotonicity arguments. Thus [the approach:] simply
    %   loop over the first successive n of these intervals,subdividing
    %   the k0-th interval appropriately. In each interval, use fzero to
    %   find the unique root there.
    %
    %   TODO
    %   - ..
    %
    %   REFERENCES
    %   [1] Predicting emissions of volatile and semivolatile organic 
    %       compounds from building materials: A review; Zhe Liu, Wei Yea,
    %       John C. Little; Building and Environment 64 (2013) 7-25.
    %   [2] http://www.mathworks.com/matlabcentral/answers/158191-
    %       positive-roots-of-trigonometric-equation
    %   [3]                       1                            1
    %       g(qn) = ------------------------------ = ---------------------
    %                beta*K + (alpha-qn.^2)*K/Bim        beta*K
    %                ----------------------------     ----------- + K/Bim
    %                         alpha-qn.^2             alpha-qn.^2
    %
    
    
    % --- 1. Initialize, check validity of parameter(s).
    %   So far, only the product K*beta is concerned. We need to have
    %   g(qn) strictly monotonic decreasing, which requires K*beta > 0,
    %   as dg/dqn = -2*K*beta*qn/()^2 (where the denominator is positive),
    %   we have sign(dg/dqn) = -sign(K*beta).
    
    internals = [] ;
    hasCorrectedErrors = false ;
    
    alpha = Q * L^2 / (D * V) ;
    beta  = A * L / V ;
    Bim   = hm * L / D ;
    
    if K * beta <= 0
        error( ['Approach based on constant monotonicity invalid for ', ...
            'negative values of product K*beta.'] ) ;
    end
    
    if nargin < 9 || isempty( correctErrors )
        correctErrors = true ;
    end
    if nargin < 10 || isempty( bisCntMax )
        bisCntMax = 1000 ;
    end
    if nargin < 11 || isempty( tol )
        tol = 1e-10 ;
    end
    
    % - Compute pole.
    %   The pole is the value of qn which nullifies the denominator
    %   of g. It is given by K*beta + (alpha-qn_pole^2)*K/Bim = 0, which
    %   implies qn_pole = +sqrt(alpha+beta*Bim) as we are interested in
    %   positive roots only.
    pole = sqrt( alpha + beta * Bim ) ;

    
    % --- 2. Compute interval boudaries, insert pole if relevant.

    % - Pre-compute interval boundaries.
    bnd = [0, pi/2 + pi * (0 : nRoots-1)] ;
    
    % - Insert pole if within boundaries and not too close from any
    %   boundary. If too close (with a few EPS) from a particular boundary,
    %   it is not possible to define a meaningful interval between them
    %   where to search for a zero and we discard the pole.
    signDiffPole = sign( bnd - pole ) ;
    if ~all( signDiffPole == signDiffPole(1) )
        ix = find( signDiffPole == 1, 1 ) ;
        bnd(ix+1:end) = bnd(ix:end-1) ;
        bnd(ix) = pole ;
    end
    
    % - Update boundaries so they are an eps aside from singularities.
    lb = bnd(1:end-1) ;
    lb = lb + 2*eps(lb) ;
    ub = bnd(2:end) ;
    ub = ub - 2*eps(ub) ;
    
    
    % --- 3. Apply a vector version of the bisection approach.
    
    bisCnt = 0 ;
    
    % - Define vectors for recording corrections.
    isCorrectedToLb = false( size( lb )) ;
    isCorrectedToUb = false( size( lb )) ;
    isCorrected = false( size( lb )) ;
    
    % - Define first mid points and errors.
    mid = (lb + ub) / 2 ;   
    err = abs( diff_fg( mid )) ;
    
    % - Define veector of logicals that flag intervals larger than EPS.
    %   The purpose is that we can stop bisecting once all intervals are
    %   equal to EPS, so we must keep track of intervals that remain
    %   larger than EPS.
    isIntervalWidthGtEps = ub - lb > eps( ub ) ;
    
    % - Iterate as long as maxIt not reached AND there are intervals whose
    %   width is greater than EPS and diff_fg(qn) greater than tolerance.
    while bisCnt <= bisCntMax && any( err(isIntervalWidthGtEps) > tol )
        bisCnt = bisCnt + 1 ;
        
        % - Compute diff_fg at interval boundaries and check that the
        %   signs differ. The SIGN function returns 0 when its arg is 0
        %   (the target). The criterion for error (due to same sign) is
        %   hence that both signs are either -1 or +1.
        diff_fg_lb = diff_fg( lb ) ;
        diff_fg_ub = diff_fg( ub ) ;
        hasSameNzSign = abs( sign(diff_fg_lb) + sign(diff_fg_ub) ) == 2 ;
        
        % - Unflag elements that were already corrected.
        hasSameNzSign = hasSameNzSign & ~isCorrected ;
        
        if any( hasSameNzSign )            
            % Given the fact that f is strictly increasing on intervals and
            % g strictly decreasing, this happens when g is extremely large
            % (with either sign) and intersects between f and g happen at
            % locations where f "explodes" and is extremly steep, within
            % real and actual/numeric interval boundaries (where actual
            % = real +/- 2*eps). In these situations diff_fg of both lower
            % and upper bounds are on the same side of the intersect and
            % have hence the same sign.
            if ~correctErrors
                %ixErrors = find( hasSameNzSign ) ;
                %for ix = ixErrors
                %   fprintf( 'Invalid condition on signs on interval %d : f([%.2g,%.2g]) = [%.2g,%.2g].\n', ...
                %       ix, lb(ix), ub(ix), diff_fg_lb(ix), diff_fg_ub(ix) ) ;
                %end
                error( 'Condition at interval(s) boundaries (opposite signs) not verified.' ) ;
            end

            % Instead of throwing an error, we define the zero at the
            % relevant interval boundary, because we know that it is
            % within 2*eps of this boundary. We also indicate that this
            % correction was made by setting hasCorrectedErrors (3rd
            % output parameter) to true.
            hasCorrectedErrors = true ;
            
            % - Define g because the test whether to correct to lb or to ub
            %   is based on the value of g. It would be cleaner to have it
            %   defined at only one place and call it from function
            %   diff_fg, but it would be slower -> ok as this clause is an
            %   error management case that may not be exploited.
            g = @(qn) (alpha - qn.^2) ./ (K * beta + (alpha - qn.^2) * K / Bim) ;
            isPos = g( lb ) > 0 ;
            
            % - Correct to ub (move lb to ub) and record operation.
            flag = hasSameNzSign & isPos ;
            lb(flag) = ub(flag) ;
            isCorrectedToUb(flag) = true ;
            
            % - Correct to lb (move ub to lb) and record operation.
            flag = hasSameNzSign & ~isPos ;
            ub(flag) = lb(flag) ;
            isCorrectedToLb(flag) = true ;
            
            % - Update record of all corrections.
            isCorrected = isCorrected | hasSameNzSign ;
            
            % - Quit bisection loop if all intervals are corrected. because
            %   none is valid for bisection anymore in this case.
            if all( isCorrected )
                break
            end
        end
        
        % - Update boundaries.
        %   NOTES
        %   - If product between diff_fg evaluated at lb and mid is
        %     strictly positive, the new interval must be [mid,ub].
        %     If it is negative or null (one of them is the zero),
        %     the new interval must be [lb,mid].
        %   - What we implement below is
        %       ub(~isMidUb) = mid(~isMidUb) ;
        %       lb(isMidUb) = mid(isMidUb) ;
        %     but in a fashion that exploits the fact that multiplication
        %     is faster than indexing (in this case).
        isMidUb = diff_fg_lb .* diff_fg( mid ) > 0 ;
        isLbMid = ~isMidUb ; 
        % Update ub where new interval is [lb, mid].
        ub = ub .* isMidUb + mid .* isLbMid ;
        % Update lb where new interval is [mid, ub].
        lb = lb .* isLbMid + mid .* isMidUb ;
        
        % - Update flags "interval width > EPS".
        isIntervalWidthGtEps = ub - lb > eps( ub ) ;
        
        % - Compute mid points and errors.
        mid = (lb + ub) / 2 ;        
        err = abs( diff_fg( mid )) ;
        
    end
    
    % - Get point with smallest error among lb, mid, and ub.
    buffer = [lb; mid; ub] ;
    [~, rowIx] = min( abs( diff_fg( buffer ))) ;
    qn = buffer( sub2ind( size(buffer), rowIx, 1:numel(mid) )) ;
    
    % - Finish with linear interpolation for low maxIt?
    %   Not relevant in this context.
    
    
    % --- 4. Build `internals` output if queried.
    
    if nargout > 2
        internals.diff_fg = @diff_fg ;
        internals.error = err ;
        internals.lb = lb ;
        internals.mid = mid ;
        internals.ub = ub ;
        internals.isIntervalWidthGtEps = isIntervalWidthGtEps ;
        internals.bisectionCnt = bisCnt ;
        internals.isCorrected = isCorrected ;
        internals.isCorrectedToLb = isCorrectedToLb ;
        internals.isCorrectedToUb = isCorrectedToUb ;
    end
    
    
    % === Nested functions.
    
    function dif = diff_fg( qn )
        dif = qn .* tan(qn) - (alpha - qn.^2) ./ (K * beta + (alpha - qn.^2) * K / Bim) ;
    end

end
