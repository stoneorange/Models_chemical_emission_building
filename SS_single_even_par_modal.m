function [Y,As,Mf_emit] = SS_single_even_par_modal(Dm, Km, dm, N, hm, Am, V, Q, Kpa, TSP, T, IC)
%% State-space (SS) discretization of a 1-D diffusional source
% Only a single diffusional source material. 
% Even discretization. 
% Ref: W Yan, Y Zhang, X Wang, 2009, 
  ...Simulation of VOC emissions from building materials by using the state-space method
  ...Building andEnvironment 44(2009)471–478
% Generates result Y, which is a (N+1)*(Nt+1) matrix that contains the...
  ... concentration at each spatial node from the bottom of source material
  ...(first row of Y) and in indoor air (last row of Y) at each time step.
  

%% Calculate parameters
dy = dm/(N-1);   % Thickness of each layer, "delta y" (m)


%% Construct the system of ODEs
P = [ [1; zeros(N-3,1)], [+diag(ones(N-3,1),-1) - 2*eye(N-2) + diag(ones(N-3,1),1)], [zeros(N-3,1); 1] ];

As = [ -2*(Dm/dy^2), 2*(Dm/dy^2), zeros(1,N-1) ;
        (Dm/dy^2) * P , zeros(N-2,1) ;
        zeros(1,N-2) , 2*(Dm/dy^2) , -2*(Dm/dy^2+hm/(Km*dy)) , 2*hm/dy/(1+Kpa*TSP) ;
        zeros(1,N-1) , Am*hm/(Km*V) , -(Q/V+Am*hm/V/(1+Kpa*TSP)) ];
Bs = [zeros(N,1); Q/V ];
Cs = eye(N+1);
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
    Mf_emit(i) = hm*Am * (q0'.*(X(end-1,:)/Km - X(end,:)/(1+Kpa*TSP))) * ((exp(Ve*T(i)) - exp(0))./Ve) / M0;
end


end

