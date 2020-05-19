 function [Y,As,Mf_emit,Mf_emit_decom,C_air_avg] = MOL_single_uneven_q_2_par_modal(Dm, Km, dm, N, hm, Am, V, Q, Kpa, TSP, T, IC)
%% Method of Lines (MOL) discretization of a 1-D diffusional source
% Only a single diffusional source material. 
% Uneven discretization with a thickness ratio of 2 between two adjacent layers.
% Apdated from Henri Gavin's treatment of boundary conditions. 
% Ref: Lei Huang, MOL-uneven-notes-LH-20190814.pdf
% Initial condition is a (N-1)*1 vector in which the first row is 0 and the
  ... other rows are C0.
% Generates result Y, which is a (N+1)*(Nt+1) matrix that contains the...
  ... concentration at each spatial node from the indoor air (first row of Y) 
  ... to the bottom of source material (last row of Y) at each time step.


%% Calculate parameters
q = 2;
dy = dm*(1-q)/(1-q^(N-1));   % Thickness of the first layer, "delta y" (m)


%% Construct the system of ODEs

Z2 =  3/2*(Dm/dy) / ( 4/3*Dm/dy + hm/Km );
Z3 = -1/6*(Dm/dy) / ( 4/3*Dm/dy + hm/Km );
Z0 = hm/(1+Kpa*TSP)   / ( 4/3*Dm/dy + hm/Km );

R2 =  2/3*Z2 - 1;
R3 =  2/3*Z3 + 1/3;
R0 =  2/3*Z0;

S2 = (Am*Dm/(dy*V)) * ( -4/3*Z2 + 3/2 );
S3 = (Am*Dm/(dy*V)) * ( -4/3*Z3 - 1/6 );
S0 = (Am*Dm/(dy*V)) * ( -4/3*Z0 ) - Q/V;


% discrete time dynamics matrix
As = zeros(N-1, N-1);
As(1,1:3) = [S0, S2, S3];
As(2,1:3) = [Dm/(dy^2)*R0, Dm/(dy^2)*R2, Dm/(dy^2)*R3];
for k = 3:N-2
    As(k, k-1:k+1) = [2*Dm/(3*2^(2*k-4)*dy^2), -3*Dm/(3*2^(2*k-4)*dy^2), Dm/(3*2^(2*k-4)*dy^2)];
end
As(N-1,N-2:N-1) = [6/5*Dm/(3*2^(2*N-6)*dy^2), -6/5*Dm/(3*2^(2*N-6)*dy^2)];

Bs = [ zeros(N-2,1); Q/V ];

Cs = [ 1, zeros(1,N-2);
       Z0, Z2, Z3, zeros(1,N-4);
       zeros(N-2,1), eye(N-2,N-2);
       zeros(1,N-3), -4/5, 9/5 ];

Ds = zeros(N+1,1);


%% Solve the system of ODEs   
[X,Ve] = eig(As);     % solve the eigen-value probem, X=eig-vecs, V=eig-vals
% for i = 1:size(As,1)
%     if Ve(i,i) > 0
%         Ve(i,i) = -1e30;
%     end
% end
Ve = diag(Ve);        % a vector of eig-vals
q0 = X \ IC;        % initial conditions in modal coordinates
xm = NaN(size(As,1),length(T));      % initialize modal superposition solution array
for j = 1:length(T)         % superposition of modal responses at each point in time
  xm(:,j) = X*diag(exp(Ve*T(j)))*q0;  
end

Y = Cs*xm;


%% Calculate Mf_emit using the integration of (Cs/Km-Ca/(1+Kpa*TSP))*hm*Am*dt
M0 = IC(end)*Am*dm;
Mf_emit = NaN(1, length(T));
for i = 1:length(T)
    Mf_emit(i) = hm*Am * (q0'.*(Z2/Km*X(2,:) + Z3/Km*X(3,:) + (Z0/Km-1/(1+Kpa*TSP))*X(1,:))) * ((exp(Ve*T(i)) - exp(0))./Ve) / M0;
end
Mf_emit_decom = NaN(size(As,1),length(T));
for i = 1:length(T)
    Mf_emit_decom(:,i) = hm*Am * (q0'.*(Z2/Km*X(2,:) + Z3/Km*X(3,:) + (Z0/Km-1/(1+Kpa*TSP))*X(1,:)))' .* ((exp(Ve*T(i)) - exp(0))./Ve) / M0;
end


%% Calculate the average air concentration over the simulation period
for i = 1:length(T)
    C_air_avg = q0' .* X(1,:) ./ Ve' * (exp(Ve*T(i))-exp(0)) / (T(i)-0) ;
end


end

