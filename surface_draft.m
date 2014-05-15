clear all
clc


a = 1;       % X Side length [m]
b = 1;       % Y Side length [m]
t = 0.1;       % Shell Thickness [m]
sigma = 1;
F = sigma*b*t;
ele_type = 'AHMAD4';
mesh = Factory.ShellMesh(ele_type,[4,2],[a,b,t]);

f_surface = @(x,y,z) near(x,a);
p_surface = find(mesh.find_nodes(f_surface));

q = [sigma 0 0 0 0]';

surface_fun = @(element,sc,sv) Physics.apply_surface_load(element,2,q,sc,sv);
L = mesh.integral_along_surface(5,0,p_surface,surface_fun)