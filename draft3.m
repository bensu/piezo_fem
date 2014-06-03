clc
clear all

%% Geometry and Mesh
a = 100e-3;             % [m] Beam's length
b = 5e-3;               % [m] Beam's width
t = 1e-3;               % [m] Beam's thickness
mesh = Factory.ShellMesh(EleType.AHMAD8,[1,2],[a,b,t]);

%% Laminate and Material
E   = 2e9;              % [Pa] Elasticity Coefficient
nu  = 0;                % Poisson coefficient
rho = 1;                % [kg/m3] Density
d13 = -0.046*(1e3);     % [(1e-3)C/m2] Piezo constant
e3  = (1e3)*0.01062e-9; % [(1e-3)C/Nm2] Permitivity

piezo_matrix = zeros(3,6);
piezo_matrix(3,1:3) = [d13 d13 0];
pzt = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);
l = 0.5;
laminate = Laminate([pzt, pzt],[l*t,(1-l)*t]); % is it /2 ???

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

%% LOADS
% All elements charged
V = 1e-3;
fem.loads.ele_vals.vals(:,1) = V;
fem.loads.ele_vals.vals(:,2) = -V;

%%

laminate.mat_num(0)
laminate.mat_num(1)
laminate.mat_num(-1)

%%

fem.solve
u_max = max(fem.dis.node_vals.vals(:,1))
v_min = min(fem.dis.node_vals.vals(:,2))
w_max = max(abs(fem.dis.node_vals.vals(:,3)))
V = fem.dis.ele_vals.vals(1,1)
