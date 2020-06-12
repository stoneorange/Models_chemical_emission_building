function [ M_source, M_air, matrix_M, C_air_avg ] = Sim_K_limited_no_sorption_par( hm,dm,Am,Km,V,Q,Kpa,TSP,v0,t )
% Simplify the single diffusion model by assuming chemical is always
% uniformly distributed in the source material.
% Such a simplified model contains two coupled ordinary ODEs which can be
% solved analytically.

% Ref of solution: https://en.wikipedia.org/wiki/Matrix_differential_equation

% Construct the constant matrix M.
matrix_M     = [-hm/(dm*Km), hm*Am/V/(1+Kpa*TSP); ...
         hm/(dm*Km), -hm*Am/V/(1+Kpa*TSP)-Q/V ];

% Get the eigenvalues and eigenvectors of the constant matrix M.
[eig_vector,eig_value] = eig(matrix_M);
eig_value = diag(eig_value);
eig_vector_1 = eig_vector(:,1);
eig_vector_2 = eig_vector(:,2);

% Calculate constants c1, c2, ..., cn for the solution
cons = eig_vector^(-1)*v0;  % v0 is the vector of initial values.

% Calculate the results vector
Results = cons(1)*eig_vector_1*exp(eig_value(1)*t)+cons(2)*eig_vector_2*exp(eig_value(2)*t);
M_source = Results(1,:);
M_air = Results(2,:);

% Calculate the average C_air over the simulation period
% C_air_avg = integrate(C_air(t)dt) / V / (t_end - t_0)
C_air_avg = (cons(1)*eig_vector_1(2)*(exp(eig_value(1)*t)-exp(0))/eig_value(1) + ...
    cons(2)*eig_vector_2(2)*(exp(eig_value(2)*t)-exp(0))/eig_value(2)) / V / (t-0) ;

end

