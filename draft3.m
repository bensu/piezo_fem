clear all
clc
%% PRELIMINARY ANALYSIS AND PARAMETERS
% Beam Theory analysis
a  = 0.5;        % X Side length [m]
b  = 0.05;       % Y Side length [m]
t  = 1e-3;       % Shell Thickness [m]
E  = 69e9;       % Elasticity Modulus [Pa]
nu = 0.3;         % Poisson Coefficient
rho = 2700;     % Density [kg/m3]
A = b*t;        % Area [m2]
I = b*(t^3)/12; % Second moment of Area [m4]
F = 5e6;
% Expected Frequencies
c = sqrt(E*I/(A*rho*a^4));
lambda = [1.875,4.694,7.885]';
expected_f = (lambda.^2)*c/(2*pi);

%% FEM and MESH
% Elements along the side
dofs_per_node = 5;
dofs_per_ele = 0;
n = 6;
mesh = Factory.ShellMesh(EleType.AHMAD8,[2*n,n],[a,b,t]);
laminate = Laminate(Material(E,nu,rho),t);
M = @(element) Physics.M_Shell(element,laminate,3);
K = @(element) Physics.K_Shell(element,laminate,3);
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
f_border = (@(x,y,z) (abs(x-a) < tol));
border = find(mesh.find_nodes(f_border));
q = [0 0 F 0 0]';
load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
    element,2,q,sc,sv);
L = mesh.integral_along_surface(dofs_per_node,dofs_per_ele, ...
    border,load_fun);
fem.loads.node_vals.dof_list_in(L);

%% Assembly
F = ~fem.bc.all_dofs;
S = fem.S;
M = fem.M;
R = fem.loads.all_dofs;

%% Damping
alpha = 0.0012;     % Return 0.0137 for damping_first_mode
beta  = 0.05; 
C = alpha*S + beta*M;

%% EigenValues
number_of_modes = 3;
n_dofs = mesh.n_dofs(dofs_per_node,dofs_per_ele);
[V,S2] = eigs(S(F,F),M(F,F),number_of_modes,'SM');
S2 = diag(S2);
M2 = diag(V'*M(F,F)*V);
for i = 1:size(V,2)
    V(:,i) = V(:,i) / (norm(V(:,i))*sqrt(M2(i)));
end
C2 = diag(V'*C(F,F)*V);
R2 = V'*R(F);

Z0 = R2 ./ S2;

%% Function handle for solver

A = [   zeros(number_of_modes)    eye(number_of_modes);
        -diag(S2)                           -diag(C2) ];
B = zeros(number_of_modes*2,1);
    
diff_z = @(t,z) (A*z + B);
z0     = zeros(number_of_modes*2,1);
[T,Z] = ode45(diff_z,[0 12],[Z0 ; zeros(number_of_modes,1)]);
plot(T,Z(:,1))

%% Modal Decomposition
dt = 5e-3;

% S2 = diag(V'*S(F,F)*V);
aux = V'*C(F,F)*V;
damping_first_mode = aux(1,1)/(2*sqrt(S2(1,1)))

%% Rebuild w(t) for one node

node_id = border(1);
node_dofs = mesh.node_dofs(dofs_per_node,node_id);
w_t = V(node_dofs(3),:)*Z(:,1:3)';
plot(T,w_t)

%% Rebuild D, final value
D = fem.compound_function(0);
aux = D.all_dofs;
aux(F) = V*Z(1,1:3)';
D.dof_list_in(aux);
dz = max(abs(D.node_vals.vals(:,3)));
%% Static Solution
fem.solve();
max_dz = max(fem.dis.node_vals.vals(:,3));
%% Compare
error = abs((dz - max_dz)/max_dz);
% testCase.verifyTrue(100*error < 1); % Error < 1%