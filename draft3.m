clc
clear all

%% Parameters
s = 1e8;        % [Pa] Applied and expected tension
E = 123e9;      % [Pa] Elastic Modulus
nu = 0.3;       % Poisson coefficient
rho = 1;        % [kg/m3] Density
t = 0.01;       % [m] Thickness
a = 0.24;       % [m] Plate length x-side
b = 0.1;        % [m] Plate length y-side
d13 = -5*(1e3);         % Piezo constant [(1e-3)C/m2]
e3 = (12.5e-9)*(1e3);   % Permitivity constant [(1e-3)C/Nm2]

%% Expected Values
% Solve constitutive equation by assuming sigma_yy = 0; D_z = 0
% Constitutive Matrix
aux = E/(1-nu^2)*[1   nu
    nu  1];
C = blkdiag(aux,-e3);
C(3,1:2) = -d13;
C(1:2,3) = -d13;
S = [s 0 0]';   % Applied tension and D_z
x = C \ S;  % Solution
% Solution decomposition
expected_u = x(1)*a;
expected_v = x(2)*b;
expected_V = x(3)*t/2;  % Not clear why is the /2 needed

%% FEM and MESH
mesh = Factory.ShellMesh(EleType.AHMAD4,[4,2],[a,b,t]);
piezo_matrix = zeros(3,6);
piezo_matrix(3,1:3) = [d13 d13 0];
material = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);
laminate = Laminate(material,t);
% Create the FemCase
dofs_per_node = 5;
dofs_per_ele = 1;
K = @(element) Physics.K_PiezoShell(element,laminate,2);
physics = Physics(dofs_per_node,dofs_per_ele,K);
fem = FemCase(mesh,physics);

%% BC
% Simply supported
tol = 1e-9;
% One edge restricted on x direction
f_edge = (@(x,y,z) (abs(x) < tol));
base = mesh.find_nodes(f_edge);
fem.bc.node_vals.set_val(base,[true false false false false]);
% One fixed point
f_corner = (@(x,y,z) (norm([x y z]) < tol));
corner = mesh.find_nodes(f_corner);
fem.bc.node_vals.set_val(corner,true);

%% LOADS
% Load the other side
x1_edge = (@(x,y,z) (abs(x-a) < tol));
border = find(mesh.find_nodes(x1_edge));
s = [s 0 0 0 0]';
load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
    element,2,s,sc,sv);
L = mesh.integral_along_surface(dofs_per_node,0, ...
    border,load_fun);
fem.loads.node_vals.dof_list_in(L);

%% Solve
fem.solve();

%%
dofs = mesh.all_eles_dofs(dofs_per_node,dofs_per_ele,1);
D = fem.dis.all_dofs;
ele_dis = D(dofs);

size(ele_dis)


u_max = max(fem.dis.node_vals.vals(:,1));
v_min = min(fem.dis.node_vals.vals(:,2));
V = fem.dis.ele_vals.vals(1,1);