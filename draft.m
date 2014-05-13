d = -5;
s = 1e8;
E = 123e9;
e = 12.5e-9;
t = 0.01;
a = 0.24;

B = [s 0 0]';
A = [   E 0 -d;
        0 E -d;
        -d -d -e];

x = A \ B;

ux = x(1)*a
V = t*x(3)/2

%% Generate Mesh 1

a = 0.24;
b = 0.12;
t = 0.01;

mesh = Factory.ShellMesh('AHMAD4',[4,2],[a,b,t]);

%% Generate Mesh 2
coords = [ 0.04 0.02    0;
           0.08 0.08    0;
           0.18 0.03    0;
           0.16 0.08    0;
           0    0       0;
           0    b       0;
           a    0       0;
           a    b       0];

connect = [ 5 1 2 6;
            5 7 3 1;
            1 3 4 2;
            6 2 4 8;
            7 8 4 3];

n_node = size(coords,1);
thickness = t*ones(1,n_node);
mesh2 = Mesh('AHMAD4',coords,connect,thickness);

%% Compare

a1 = mesh.nodes_per_ele
a2 = mesh2.nodes_per_ele

