function [Y,As,Mf_emit] = MOL_single_even_par_modal(Dm, Km, dm, N, hm, Am, V, Q, Kpa, TSP, T, IC)
%% Method of Lines (MOL) discretization of a 1-D diffusional source
% A single diffusional source material. 
% Even discretization. 
% Huang's even method.
% Initial condition is a (N-1)*1 vector in which the last row is 0 and the
  ... other rows are C0.
% Generates result Y, which is a (N+1)*(Nt+1) matrix that contains the...
  ... concentration at each spatial node from the bottom of source material
  ...(first row of Y) and in indoor air (last row of Y) at each time step.


%% Calculate parameters
dy = dm/(N-1);   % Thickness of each layer (m)
L = Am/V;


%% Construct the system of ODEs
P = +diag(ones(N-3,1),-1) - 2*eye(N-2) + diag(ones(N-3,1),1);

Z3 =  -0.5*(Dm/dy) / ( 1.5*Dm/dy + hm/Km );
Z2 =   2.0*(Dm/dy) / ( 1.5*Dm/dy + hm/Km );
Z0 =  hm/(1+Kpa*TSP) / ( 1.5*Dm/dy + hm/Km );

R3 =  Z3 + 1.0;
R2 =  Z2 - 2.0;
R0 =  Z0;

S3 = -(L*Dm/dy) * ( 1.5*Z3 + 0.5 );
S2 = -(L*Dm/dy) * ( 1.5*Z2 - 2.0 );
S0 = -(L*Dm/dy) * ( 1.5*Z0 ) - Q/V;

P(1,:) = (2/3) * [ -1 , 1 , zeros(1,N-4) ];  % Neumann BC at y=0, bottom of the source material.

P(N-2,:) = [ zeros(1,N-4), R3, R2 ];  % Dmirichmlet BC at y=dm, surface of the source material. 

% discrete time dynamics matrix
As = [ (Dm/dy^2) * P , [ zeros(N-3,1) ; Dm/dy^2*R0 ] ;
        zeros(1,N-4) , S3 , S2 , S0 ] ;
Bs = [ zeros(N-2,1); Q/V ];
Cs = [ 4/3 , -1/3 , zeros(1,N-3) ;
      [ eye(N-2,N-1) ;
       zeros(1,N-4) , Z3 , Z2 , Z0 ;
       zeros(1,N-2) , 1 ] ];
Ds = zeros(N+1,1);


%% Solve the system of ODEs   
[X,Ve] = eig(As);     % solve the eigen-value probem, X=eig-vecs, V=eig-vals
Ve = diag(Ve);        % a vector of eig-vals
q0 = X \ IC;        % initial conditions in modal coordinates
xm = NaN(size(As,1),length(T));      % initialize modal superposition solution array
for j = 1:length(T)         % superposition of modal responses at each point in time
  xm(:,j) = X*diag(exp(Ve*T(j)))*q0;  
end

Y = Cs*xm;


%% Calculate Mf_emit using the integration of (Cs/Km-Ca/(1+Kpa*TSP))*hm*Am*dt
M0 = IC(1)*Am*dm;
Mf_emit = NaN(1, length(T));
for i = 1:length(T)
    Mf_emit(i) = hm*Am * (q0'.*(Z2/Km*X(end-1,:) + Z3/Km*X(end-2,:) + (Z0/Km-1/(1+Kpa*TSP))*X(end,:))) * ((exp(Ve*T(i)) - exp(0))./Ve) / M0;
end


end

