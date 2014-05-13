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

x = As \ Bs;

ux = x(1)*a
V = t*x(3)/2