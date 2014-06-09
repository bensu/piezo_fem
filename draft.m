clear all
clc

%% Represent a spring-mass system using shells

%% Parameters
a = 1;
L = 2;
t = 1;
E = 1;       % Elasticity Modulus [Pa]
nu = 0;      % Poisson Coefficient
rho = 0;     % Density [kg/m3]
A = a*t;     % Area [m2]
mass = 1;
%% Expected Values
k = E*A/L;
expected_omega = sqrt(k/mass);

%% FEM & Mesh
dofs_per_node = 5;
dofs_per_ele = 0;
n = 6;
mesh = Factory.ShellMesh(EleType.AHMAD4,[1,1],[a,L,t]);

% Add point mass
tol = 1e-9;
f_edge = @(x,y,z) (abs(y - L) < tol);
edge_nodes = find(mesh.find_nodes(f_edge));
mass_values = ones(size(edge_nodes));
mass_values = mass*mass_values/sum(mass_values);
mesh = mesh.add_point_mass(edge_nodes,mass_values);

laminate = Laminate(Material(E,nu,rho),t);
M = @(element) Physics.M_Shell(element,laminate,3);
K = @(element) Physics.K_Shell(element,laminate,3);
physics = Physics.Dynamic(dofs_per_node,dofs_per_ele,K,M);
fem = FemCase(mesh,physics);

%% Boundary Conditions
% Restrict all non-y movement
fem.bc.node_vals.vals(:,[1 3:end]) = true;
f_clamped = @(x,y,z) (abs(y) < tol);
clamped_noes = mesh.find_nodes(f_clamped);
fem.bc.node_vals.set_val(clamped_noes,true);

%% Dynamic Problem

mode_number = 2;
[Z,D] = fem.eigen_values(mode_number);
omega = sqrt(diag(D))

omega(1) - expected_omega(1)

%%

Z(1).node_vals.vals(:,2)