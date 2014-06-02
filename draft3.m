

m1 = Material(1,1,1);
m2 = Material(2,1,1);
m3 = Material(3,1,1);

laminate    = [m1,m2,m3];
thickness   = [1,3,4];


t_sum   = sum(thickness);

zeta_vals = interp1([0 t_sum/2 t_sum],[-1 0 1],cumsum(thickness) ...
                                        ,'linear','extrap')

zeta = 0;

m = laminate(find(zeta <= zeta_vals,1))