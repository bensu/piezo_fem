clc
clear all
type = EleType.H8;
n_elements = [1 1 1];
sides = [1 1 1];
mesh = Factory.BrickMesh(type,n_elements,sides);
e = mesh.ele(1);
all(0.5*ones(3,1) == diag(e.jacobian(0,0,0)));


