clear all
clc
%% PRELIMINARY ANALYSIS AND PARAMETERS
% Beam Theory analysis
a = 0.5;    % X Side length [m]
b = 0.05;   % Y Side length [m]
t = 1e-3;   % Shell Thickness [m]
E = 69e9;	% Elasticity Modulus [Pa]
nu = 0.33;   % Poisson Coefficient
rho = 2700;	% Density [kg/m3]
A = b*t;        % Area [m2]
I = b*(t^3)/12; % Second moment of Area [m4]

% Expected Frequencies
c = sqrt(E*I/(A*rho*a^4));
lambda = [1.875,4.694,7.885]';
beam_freq = @(L) ((lambda.^2)*sqrt(E*I/(A*rho*L^4))/(2*pi));
expected_f = beam_freq(a)
expected_f_400 = beam_freq(0.4);   % For Project Report
expected_f_600 = beam_freq(0.6);

%% FEM and MESH
% Elements along the side
dofs_per_node = 5;
dofs_per_ele = 0;
n = 10;
m = 20;
lam = Laminate(Material(E,nu,rho),t);
mesh = Factory.ShellMesh(EleType.AHMAD8,lam,[m,n],[a,b,t]);

M = @(element) Physics.M_Shell(element,3);
K = @(element) Physics.K_Shell(element,3);
physics = Physics.Dynamic(dofs_per_node,dofs_per_ele,K,M);
fem = FemCase(mesh,physics);

%% BC
% Fixed End
tol = 1e-9;
x0_edge = (@(x,y,z) (abs(x) < tol));
base = mesh.find_nodes(x0_edge);
fem.bc.node_vals.set_val(base,true);

%% Calculate Frequencies
number_of_modes = 3;    % To be calculated
[Z, W] = fem.eigen_values(number_of_modes);
found_f = sqrt(diag(W))/(2*pi)
error = abs(found_f(1) - expected_f(1)) / expected_f(1);

%%  Plot modal shapes
f_side = @(x,y,z) (y < tol);
side_nodes = mesh.find_nodes(f_side);
x_coords = mesh.coords(side_nodes);
for n = 1:number_of_modes
    subplot(number_of_modes,1,n);
    hold on
    % Normalize
    z = Z(n).node_vals.vals(side_nodes,3);
    [val, loc] = max(abs(z));
    title(['Modal Shape ',num2str(n)])
    xlabel('x [m]')
    ylabel('Normalized Amplitude [m]')
    plot(x_coords,sign(z(loc))*z/val,'ko')
    plot(x_coords,beam_shapes(a,n,x_coords))
    legend('Normalized EigenVector','Beam Shape Function');
    hold off
end

%% Antinode for mode 2
bs = beam_shapes(a,2,x_coords);
index = find(x_coords>0.2);
[min_bs,l] = min(abs(bs(index)));
hold on
plot(x_coords,bs)
plot(x_coords(index(l)),bs(index(l)),'o')
hold off
