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
d13 = -5*(1e3);         % Piezo constant [(1e-3)C/m2]
e3 = (12.5e-9)*(1e3);   % Permitivity constant [(1e-3)C/Nm2]
% Expected Frequencies
c = sqrt(E*I/(A*rho*a^4));
lambda = [1.875,4.694,7.885]';
expected_f = (lambda.^2)*c/(2*pi);

%% FEM and MESH
% Elements along the side
dofs_per_node = 5;
dofs_per_ele = 1;
n = 6;
mesh = Factory.ShellMesh('AHMAD8',[2*n,n],[a,b,t]);
piezo_matrix = zeros(3,6);
piezo_matrix(3,1:3) = [d13 d13 0];
material = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);
M = @(element) Physics.M_Shell(element,material,3);
K = @(element) Physics.K_PiezoShell(element,material,2);
physics = Physics.Dynamic(dofs_per_node,dofs_per_ele,K,M);
fem = FemCase(mesh,physics);

%% BC
% Fixed End
tol = 1e-9;
x0_edge = (@(x,y,z) (abs(x) < tol));
base = mesh.find_nodes(x0_edge);
fem.bc.node_vals.set_val(base,true);
% All grounded but one
grounded = true(mesh.n_ele,1);
grounded(5) = false;    % All except element 5
fem.bc.ele_vals.set_val(grounded,true);
fem.bc.ele_vals.vals

%% LOADS
% Load the other side
f_border = (@(x,y,z) (abs(x-a) < tol));
border = find(mesh.find_nodes(f_border));
q = [0 0 F 0 0]';
load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
    element,2,q,sc,sv);
L = mesh.integral_along_surface(dofs_per_node,0,border,load_fun);
fem.loads.node_vals.dof_list_in(L);

%% Assembly
F = ~fem.bc.all_dofs;
S = fem.S;
M = fem.M;

mech_dofs = mesh.master_node_dofs(dofs_per_node,1:mesh.n_nodes,1:dofs_per_node);
elec_dofs = mesh.master_ele_dofs(dofs_per_node,dofs_per_ele, ...
                                            1:mesh.n_ele,dofs_per_ele);
                                        
Muu = M(F(mech_dofs),F(mech_dofs));
Kuu = S(F(mech_dofs),F(mech_dofs));
Kpp = S(F(elec_dofs),F(elec_dofs));
Kup = S(F(mech_dofs),F(elec_dofs));


%% EigenValues
n_dofs = mesh.n_dofs(dofs_per_node,dofs_per_ele);
n_modes = 2;
[V,D] = eigs(Kuu,Muu,n_modes,'SM');
M2 = diag(V'*Muu*V);
for i = 1:size(V,2)
    V(:,i) = V(:,i) / (norm(V(:,i))*sqrt(M2(i)));
end

%% Modal Decomposition
dt = 5e-3;
R = fem.loads.all_dofs;
C = sparse(eye(length(mech_dofs))/10000);
S2 = diag(V'*Kuu*V);
C2 = diag(V'*C(F(mech_dofs),F(mech_dofs))*V)/(2*dt);
R2 = V'*R(F(mech_dofs));

%% Time integration
steps = 300;
Z = zeros(size(V,2),steps+1);
K2T = S2 - 2/(dt^2);
MC1 = 1./(dt^-2 + C2);
MC2 = dt^-2 - C2;

%% Control

P = V'*Kup;

A = [   zeros(size(diag(S2))) eye(size(diag(C2)))
        -diag(S2)             -diag(C2)         ]
B = [ zeros(size(P));   P]

rank(ctrb(A,B))

%% Time integration

f = @(z) (
tt = 0:dt:(steps*dt);
for n = 2:steps
    Z(:,n+1) = MC1.*(R2 - K2T.*Z(:,n) - MC2.*Z(:,n-1));
end
%% Rebuild D, final value
D = fem.compound_function(0);
aux = D.all_dofs;
aux(F(mech_dofs)) = V*(Z(:,end));
D.dof_list_in(aux);
dz = max(D.node_vals.vals(:,3));
%% Static Solution
fem.solve();
max_dz = max(fem.dis.node_vals.vals(:,3));
%% Compare
error = abs((dz - max_dz)/max_dz);