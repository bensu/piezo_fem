clc
clear all

%% PRELIMINARY ANALYSIS AND PARAMETERS
% Beam Theory analysis
a = 0.5;        % X Side length [m]
b = 0.05;       % Y Side length [m]
t = 1e-3;       % Shell Thickness [m]
E = 69e9;       % Elasticity Modulus [Pa]
nu = 0;         % Poisson Coefficient
rho = 2700;     % Density [kg/m3]
A = b*t;        % Area [m2]
I = b*(t^3)/12; % Second moment of Area [m4]

F = 1;

% Expected Frequencies
c = sqrt(E*I/(A*rho*a^4));
lambda = [1.875,4.694,7.885]';
expected_f = (lambda.^2)*c/(2*pi);

%% FEM and MESH
% Elements along the side
dofs_per_node = 5;
dofs_per_ele = 0;
mesh = Factory.ShellMesh('AHMAD8',[3,6],[a,b,t]);
material = Material(E,nu,rho);
M = @(element) Physics.M_Shell(element,material,3);
K = @(element) Physics.K_Shell(element,material,3);
physics = Physics.Dynamic(dofs_per_node,dofs_per_ele,K,M);
fem = FemCase(mesh,physics);

%% BC
% Fixed End
tol = 1e-9;
x0_edge = (@(x,y,z) (abs(x) < tol));
base = mesh.find_nodes(x0_edge);
fem.bc.node_vals.set_val(base,true);

%% LOADS
% Load the other side
x1_edge = (@(x,y,z) (abs(x-a) < tol));
border = find(mesh.find_nodes(x1_edge));
q = [0 0 F 0 0]';
load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
    element,2,q,sc,sv);
L = mesh.integral_along_surface(dofs_per_node,dofs_per_ele, ...
    border,load_fun);
fem.loads.node_vals.dof_list_in(L);

%% Static Solution

fem.solve();
[max_dz, max_node] = max(fem.dis.node_vals.vals(:,3));

%% Dynamic Solution

[V,D] = fem.eigen_values(3);

V
D

%% Damping
n_dofs = mesh.n_dofs(dofs_per_node,dofs_per_ele);

F = ~fem.bc.all_dofs;
R = fem.loads.all_dofs;
R = R(F);
M = fem.M; M = M(F,F);
S = fem.S; S = S(F,F);
C = sparse(eye(n_dofs));
C = C(F,F);



%% Force

% steps = 1000000;
% dt = a*sqrt(rho/E)/1000;   % /N added for safety.
% D = zeros(size(S,1),steps+1);
% 
% M2 = M/(dt^2);
% C2 = C/(2*dt);
% MC = M2 + C2;
% 
% for n = 2:steps
%     D(:,n+1) = MC \ (R -S*D(:,n) + M2*(2*D(:,n)-D(:,n-1)) + C2*D(:,n-1));
% end

%% Compare final value
% 
% max_dz
% max(D(:,end))
% 
% %% 
% 
% t = linspace(0,dt*steps,steps+1);
% hold on
% plot(t,D((max_node+1)*dofs_per_node+3,:))
% % plot(t,max_dz)
% hold off


