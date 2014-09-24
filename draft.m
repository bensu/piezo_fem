
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

%% MESH
% Elements along the side
laminate = Laminate(Material(E,nu,rho),t);
mesh = Factory.ShellMesh(EleType.AHMAD8,laminate,[10,5],[a,b,t]);

total_mass = 0.003;
tol = 1e-9;
between = @(a,b,x) ((x >= a) && (x <= b));
f_acc = @(x,y,z) (between(0.38,0.45,x) && between(0.02,0.03,y));
accelerometer_nodes = find(mesh.find_nodes(f_acc));
total_nodes = length(accelerometer_nodes);
mass_values = total_mass*ones(total_nodes,1)/total_nodes;
mesh = mesh.add_point_mass(accelerometer_nodes,mass_values);
%% FEM
dofs_per_node = 5;
dofs_per_ele = 0;
M = @(element) Physics.M_Shell(element,3);
K = @(element) Physics.K_Shell(element,3);
physics = Physics.Dynamic(dofs_per_node,dofs_per_ele,K,M);
fem = FemCase(mesh,physics);

%% BC
% Fixed End
x0_edge = (@(x,y,z) (abs(x) < tol));
base = mesh.find_nodes(x0_edge);
fem.bc.node_vals.set_val(base,true);

mesh.plot()

%% Calculate Frequencies
[~, W] = fem.eigen_values(3);
found_f = sqrt(diag(W))/(2*pi)
error = abs(found_f(1) - expected_f(1)) / expected_f(1);
% Error should be less than 5 percent
% testCase.verifyTrue(100*error < 5);