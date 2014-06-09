clear all
clc
close all

%% Scordelis-Lo Roof Problem
% It consists of a cylindrical roof fixed on its curved ends and loaded by 
% a constant vertical load, characterizing the roof's own weight.

%% IMPORT MESH
% The meshes are imported from the /scordelis folder
n = 10;                         % Number of elements along each side, [3,5,10,20]
ele_type = EleType.AHMAD9;      % Element Type
nodes_per_ele = ele_type.n_nodes;
coords_in = load(['tests/scordelis/nodos' num2str(n) '_' ...
                          num2str(nodes_per_ele) '.txt']);
connect_in = load(['tests/scordelis/elementos' num2str(n) '_' ...
                          num2str(nodes_per_ele) '.txt']);

%% Format Input Mesh:
% Remove Extra Columns
coords_in(:,[1 end]) = [];
connect_in(:,1) = [];
% Clean and reformat connect matrix
connect_in(connect_in == 0) = [];
connect_in = reshape(connect_in,[],nodes_per_ele);

to = 0.25;  % [m] Thickness
t = to*ones(1,size(coords_in,1));
mesh = Mesh(ele_type,coords_in,connect_in,t);

%% MATERIAL
E = 4.32E8;     % [Pa] Elasticity coefficient
nu = 0.0;       % Poisson coefficient
laminate = Laminate(Material(E,nu,1),to);

%% FEM
dofs_per_node = 5;
dofs_per_ele = 0;
fun_in = @(element) Physics.K_Shell(element,laminate,3);
physics = Physics(dofs_per_node,dofs_per_ele,fun_in);
fem = FemCase(mesh,physics);

%% LOADS
q = [0 0 -360 0 0]';    % [N/m3] Applied as a volume load
fun_in = @(element) Physics.apply_volume_load(element,3,q);
R = mesh.assembly_vector(dofs_per_node,dofs_per_ele,fun_in);
fem.loads.node_vals.dof_list_in(R);

%% BC
bc = false(mesh.n_nodes,dofs_per_node);
tol = 1e-3;
% Supported edge at x = 0
f_edge_1 = @(x,y,z) (x < tol);
bc(mesh.find_nodes(f_edge_1),[2 3]) = true;
% Symmetry at y = 0
f_center = @(x,y,z) (y < tol);
bc(mesh.find_nodes(f_center),[2 4]) = true;
% Symmetry at x = 25
f_edge_2 = @(x,y,z) (x > 25 - tol);
bc(mesh.find_nodes(f_edge_2),[1 5]) = true;
fem.bc.node_vals.vals = bc;

%% SOLVE

fem.solve();
U = fem.dis.node_vals.vals;

max(abs(U(:,3))),

%% PLOT 

scale = 1;
PlotMesh(mesh.coords + scale*fem.dis.node_vals.vals(:,1:3), ...
    mesh.connect, ...
    ~fem.bc.node_vals.vals, ...
                fem.reactions.node_vals.vals);

%% Expected Solutions

% Draw_Placas(mesh.connect,coords + U(:,1:3),'b')

% para t = 2.5000 Dz = -0.0266792 (ADINA) 1/10
% -0.026570686677499 AHMAD9
% para t = 0.2500 Dz = -0.3024  (REF)   -0.302397 (ADINA) 1/100
% -0.299240246836452 AHMAD9
% para t = 0.0250 Dz = -3.21412 (ADINA) 1/1000
% -2.072900027933127 AHMAD9
% para t = 0.0025 Dz = -3.25761E+01  (ADINA) 1/10000
% -2.800824498499196 AHMAD9

