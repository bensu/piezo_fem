clc
clear all

%% Geometry and Mesh
a = 100e-3;             % [m] Beam's length
b = 5e-3;               % [m] Beam's width
t = 1;                  % [m] Beam's thickness
mesh = Factory.ShellMesh(EleType.AHMAD8,[4,2],[a,b,t]);

%% Laminate and Material
E = 2e9;                % [Pa] Elasticity Coefficient
nu = 0;                 % Poisson coefficient
rho = 1;                % [kg/m3] Density
metal = Material(E,nu,rho);
l = 0;
% laminate = Laminate([metal, metal],[l*t,(1-l)*t]); % is it /2 ???
laminate = Laminate(metal,t);

%% Physics and FEM
dofs_per_node = 5;
dofs_per_ele = 0;
K = @(element) Physics.K_Shell(element,laminate,2);
physics = Physics(dofs_per_node,dofs_per_ele,K);
fem = FemCase(mesh,physics);

%%
zeta_m = laminate.zeta_m
zeta_t = laminate.zeta_t