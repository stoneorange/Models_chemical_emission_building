function [Y,As,Mf_emit] = SS_single_uneven_Guo_par_modal(Dm, Km, dm, N, q, hm, Am, V, Q, Kpa, TSP, T, IC)
%% State-space (SS) discretization of a 1-D diffusional source
% A single diffusional source material. 
% First layer is set to alpha=1e-7 m by default and the first node is at the surface, then the material thickness (dm-1e-7) is unevenly
% discretized into (N-1) nodes situated at the middle of the (N-1) layers.
% Guo's uneven method.

% Ref: 
%  Z Guo, 2013,
%   A Framework for Modelling Non-Steady-State Concentrations of Semivolatile Organic Compounds 
%   Indoors – I: Emissions from Diffusional Sources and Sorption by Interior Surfacs
%   Indoor Built Environment 2013(22): 685-700.
%  W Yan, Y Zhang, X Wang, 2009, 
%   Simulation of VOC emissions from building materials by using the state-space method
%   Building and Environment 44(2009)471–478

%% Calculate parameters
alpha=1e-7; % Thickness of the first layer (m)
dy = (dm-alpha)*(1-q)/(1-q^(N-2));   % Thickness of the second layer (m)

%% Construct the system of ODEs 
P = zeros(N-4,N-2);
for i = 3:(N-2)
    P(i-2,i-2) = Dm/(1/2*(dy*q^(i-2))^2*(1+1/q));
    P(i-2,i-1) = -Dm*(1/(1/2*(dy*q^(i-2))^2*(1+1/q)) + 1/(1/2*(dy*q^(i-2))^2*(1+q)));
    P(i-2,i) = Dm/(1/2*(dy*q^(i-2))^2*(1+q));
end

R0 = -(Am*hm/V/(1+Kpa*TSP)+Q/V);
R1 = Am*hm/(V*Km);

S0 = hm/alpha/(1+Kpa*TSP);
S1 = -(hm/(alpha*Km) + Dm/(alpha*(alpha+1/2*dy)));
S2 = Dm/(alpha*(alpha+1/2*dy));

Z1 = Dm/(dy*(alpha+1/2*dy));
Z2 = -Dm*(1/(dy*(alpha+1/2*dy)) + 1/(dy*(1/2*dy+1/2*dy*q)));
Z3 = Dm/(dy*(1/2*dy+1/2*dy*q));

U_N2 = Dm/(1/2*(dy*q^(N-3))^2*(1+1/q));
U_N1 = -Dm*(1/(1/2*(dy*q^(N-3))^2*(1+1/q)) + 1/(1/2*(dy*q^(N-3))^2));
U_N = Dm/(1/2*(dy*q^(N-3))^2);

As = [ R0, R1, zeros(1,N-1);
       S0, S1, S2, zeros(1,N-2);
       0, Z1, Z2, Z3, zeros(1,N-3);
       zeros(N-4,2), P, zeros(N-4,1);
       zeros(1,N-2), U_N2, U_N1, U_N;
       zeros(1,N-1), Dm/(1/2*dy*q^(N-3))^2, -Dm/(1/2*dy*q^(N-3))^2 ];
Bs = [ Q/V; zeros(N,1) ];
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
M0 = IC(end)*Am*dm;
Mf_emit = NaN(1, length(T));
for i = 1:length(T)
    Mf_emit(i) = hm*Am * (q0'.*(X(2,:)/Km - X(1,:)/(1+Kpa*TSP))) * ((exp(Ve*T(i)) - exp(0))./Ve) / M0;
end


end

