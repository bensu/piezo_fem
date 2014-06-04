%% Parameters
a = 100e-3;             % [m] Beam's length
b = 5e-3;               % [m] Beam's width
t = 1e-3;               % [m] Beam's thickness
E   = 2e9;              % [Pa] Elasticity Coefficient
nu  = 0;                % Poisson coefficient
rho = 1;                % [kg/m3] Density
d13 = -0.046*(1e3);     % [(1e-3)C/m2] Piezo constant
e3  = (1e3)*0.01062e-9; % [(1e-3)C/Nm2] Permitivity
V = 1e-3;
%% Expected Values
% [Development of a Light-Weight Robot End-Effector using
% Polymeric Piezoelectric Bimorph, Tzou, 1989][Eq 20]
expected_w = -3*d13*V*(a^2)/(2*E*t^2);
%% Geometry and Mesh
mesh = Factory.ShellMesh(EleType.AHMAD8,[1,1],[a,b,t]);
%% Laminate and Material
piezo_matrix = zeros(3,6);
piezo_matrix(3,1:3) = [d13 d13 0];
pzt = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);
l = 0.5;
laminate = Laminate([pzt, pzt],[l*t,(1-l)*t]);
%% Physics and FEM
dofs_per_node = 5;
dofs_per_ele  = laminate.n_layers;
gauss_order   = 2;
K = @(element) Physics.K_PiezoShell(element,laminate,gauss_order);
physics = Physics(dofs_per_node,dofs_per_ele,K);
fem = FemCase(mesh,physics);
%% BC
tol = 1e-9;
% Clamped
f_edge = (@(x,y,z) (abs(x) < tol));
base = mesh.find_nodes(f_edge);
fem.bc.node_vals.set_val(base,true);
%% Loads - Initial Displacements
% All elements charged, initial condition
fem.initial_dis.ele_vals.vals(:,1) = V/2;
fem.initial_dis.ele_vals.vals(:,2) = -V/2;
%% Solve
fem.solve

%% Check Stresses and Strains (for Report)
ele_id = 1;
ele = mesh.ele(ele_id);
dofs = mesh.all_eles_dofs(dofs_per_node,dofs_per_ele,ele_id);
D = fem.dis.all_dofs;
ele_dis = D(dofs);
B = Physics.B_PiezoShell(ele,laminate.n_layers,1,0,0,-1);
C = laminate.PiezoMatrix(0);
format long
C*B*ele_dis
z = t/2;
sigma_x_z = d13*V/(t^2)*z;

%% Verify
expected_w;
w_max = max(abs(fem.dis.node_vals.vals(:,3)));
