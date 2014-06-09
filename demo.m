%% Bimorph Piezoelectric Problem

% A beam consisting of two piezoelectric layers is polarized, each layer in
% a different direction. This results in opposite tensions, that apply a
% moment and bend the beam.

%% Parameters
a = 100e-3;             % [m] Beam's length
b = 5e-3;               % [m] Beam's width
t = 1e-3;               % [m] Beam's thickness
E   = 2e9;              % [Pa] Elasticity Coefficient
nu  = 0;                % Poisson coefficient
rho = 1;                % [kg/m3] Density
d13 = -0.046*(1e3);     % [(1e-3)C/m2] Piezo constant
e3  = (1e3)*0.01062e-9; % [(1e-3)C/Nm2] Permitivity
V = 1e-3;               % [(1e3)V] Applied voltage
%% Expected Values
% [Development of a Light-Weight Robot End-Effector using
% Polymeric Piezoelectric Bimorph, Tzou, 1989][Eq 20]
expected_w = -3*d13*V*(a^2)/(2*E*t^2);

%% Geometry and Mesh
% Define a Mesh with element type AHMAD8 and dimensions a x b x t
n = 4;      % Number of elements along x
m = 2*n;    % Number of elements along y
mesh = Factory.ShellMesh(EleType.AHMAD8,[m,n],[a,b,t]);

%% Laminate and Material
% Define a Piezoelectric material, and a two layer laminate
piezo_matrix = zeros(3,6);
piezo_matrix(3,1:3) = [d13 d13 0];
pzt = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);
laminate = Laminate([pzt, pzt],[t/2,t/2]);
%% Physics and FEM
dofs_per_node = 5;
dofs_per_ele  = laminate.n_layers; % 2
gauss_order   = 2;
K = @(element) Physics.K_PiezoShell(element,laminate,gauss_order);
physics = Physics(dofs_per_node,dofs_per_ele,K);
fem = FemCase(mesh,physics);
%% BC
tol = 1e-9;
% Clamped
f_edge = (@(x,y,z) (abs(x) < tol));     % function returns true for nodes in the edge
base = mesh.find_nodes(f_edge);         % Nodes in the edge
fem.bc.node_vals.set_val(base,true);    % All dofs are supported -> true
%% Loads - Initial Displacements
% All elements charged, initial condition
fem.initial_dis.ele_vals.vals(:,1) = V/2;  
fem.initial_dis.ele_vals.vals(:,2) = -V/2;
%% Solve
fem.solve
%% Verify
expected_w
w_max = max(abs(fem.dis.node_vals.vals(:,3)))
%%

f_side = @(x,y,z) (abs(y)<tol);
side_nodes = mesh.find_nodes(f_side);
x = mesh.coords(side_nodes,1);
w = fem.dis.node_vals.vals(side_nodes,3);


scale = 10000;
fem.plot(scale)