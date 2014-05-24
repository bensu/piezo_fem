clc
clear all

%% PRELIMINARY ANALYSIS AND PARAMETERS
% In-plane patch test for one element
% Exactly like 5.1.8 Lentzen Thesis
a = 0.24;       % X Side length [m]
b = 0.12;       % Y Side length [m]
t = 0.01;       % Shell Thickness [m]

E	= 123e9;    % Elasticity Modulus [N/m2]
nu	= 0;        % Poisson Coefficient / Beam Theory requires nu = 0
rho = 1;        % Density [kg/m3]
e13 = -5*(1e3);       % Piezo constant [(1e-3)C/m2]
e23 = 0*(1e3);        % Piezo constant [(1e-3)C/m2], not used for this problem
e3  = 12.5e-9*(1e3);  % Permitivity constant [(1e-3)C/Nm2]

M = 2e3;      % Sum of Moments applied to points [N]
area = b*t;
M_dist = M/(area);
I = b*(t^3)/12;

% Expected Displacemente along the x and y axis
expected_dz = M*(a^2)/(2*E*I) % [m]
expected_phi = -M*a/(E*I)
expected_V  = 1.686e3   % [V]

dofs_per_node = 5;
dofs_per_ele  = 1;

%% MESH

mesh = Factory.ShellMesh('AHMAD8',[1,1],[a,b,t]);

%% MATERIAL
piezo_matrix = zeros(3,6);
piezo_matrix(3,1:3) = [e13 e13 0];
material = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);

%% FEM and MESH
% Create the FemCase
K = @(element) Physics.K_PiezoShell(element,material,2);
physics = Physics(dofs_per_node,dofs_per_ele,K);
fem = FemCase(mesh,physics);

%% BC
% Clamped since there is no poisson effect
tol = 1e-9;
x0_edge = (@(x,y,z) (abs(x) < tol));
base = mesh.find_nodes(x0_edge);
fem.bc.node_vals.set_val(base,true);
% All bottom layer is grounded
%             fem.bc.ele_vals.vals(:,1) = true;

%% LOADS
% Load the other side
x1_edge = (@(x,y,z) (abs(x-a) < tol));
border = find(mesh.find_nodes(x1_edge));
q = [0 0 0 0 -M_dist]';
load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
    element,2,q,sc,sv);
L = mesh.integral_along_surface(dofs_per_node,0, ...
    border,load_fun);
fem.loads.node_vals.dof_list_in(L);

%% SOLUTION
fem.solve();    % SIDE EFFECTS
%% TEST
% Check if the obtained values are the expected ones
fem.loads.node_vals.vals;
fem.reactions.node_vals.vals;

dz = max(fem.dis.node_vals.vals(:,3))
phi = min(fem.dis.node_vals.vals(:,5))
V   = (fem.dis.ele_vals.all_dofs)*(1e3)
%             testCase.verifyEqual(true,near(expected_dz,dz));
%             testCase.verifyEqual(true,near(expected_V,max_V));

%% PLOT
%             PlotMesh(mesh.coords + 1*fem.dis.node_vals.vals(:,1:3), ...
%                 mesh.connect, ...
%                 ~fem.bc.node_vals.vals, ...
%                 fem.reactions.node_vals.vals);