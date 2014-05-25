clc
clear all

%% PRELIMINARY ANALYSIS AND PARAMETERS
% Beam Theory analysis
a = 0.5;    % X Side length [m]
b = 0.05;   % Y Side length [m]
t = 1e-3;   % Shell Thickness [m]
E = 69e9;	% Elasticity Modulus [Pa]
nu = 0;   % Poisson Coefficient
rho = 2700;	% Density [kg/m3]
A = b*t;        % Area [m2]
I = b*(t^3)/12; % Second moment of Area [m4]

% Expected Frequencies
c = sqrt(E*I/(A*rho*a^4));
lambda = [1.875,4.694,7.885]';
expected_f = (lambda.^2)*c/(2*pi);

%% FEM and MESH
% Elements along the side
dofs_per_node = 5;
dofs_per_ele = 0;
mesh = Factory.ShellMesh('AHMAD8',[10,5],[a,b,t]);
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

%% Calculate Frequencies
[V D] = fem.eigen_values(3);
found_f = sqrt(diag(D))/(2*pi);
(near(expected_f(1)/1000,found_f(1)/1000))
