clc
clear all

%% PRELIMINARY ANALYSIS AND PARAMETERS
% Beam Theory analysis
a = 0.5;    % X Side length [m]
b = 0.05;   % Y Side length [m]
t = 1e-3;   % Shell Thickness [m]
E = 69e9;	% Elasticity Modulus [Pa]
nu = 0.3;   % Poisson Coefficient
rho = 2700;	% Density [kg/m3]
A = b*t;        % Area [m2]
I = b*(t^3)/12; % Second moment of Area [m4]

% Expected Frequencies
c = sqrt(E*I/(A*rho*a^4));
lambda = [1.875,4.694,7.885]';
f = (lambda.^2)*c/(2*pi);

%% FEM and MESH
% Elements along the side
dofs_per_node = 5;
dofs_per_ele = 0;
mesh = Factory.ShellMesh('AHMAD4',[5,2],[a,b,t]);

%%
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
q = [1 0 0 0 0]';
load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
    element,2,q,sc,sv);
L = mesh.integral_along_surface(dofs_per_node,dofs_per_ele, ...
    border,load_fun);
fem.loads.node_vals.dof_list_in(L);
fem.solve()

%% Plot

fem.plot()

%% Obtained Mass

[V D] = fem.eigen_values();
diag(D);

