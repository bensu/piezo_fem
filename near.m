function bool_out = near(x,y)
% bool_out = near(x,y)
% Takes two real values and checks if they are equal according to an 
% absolute tolerance (Atol) and a relative tolerance (Rtol)
% Should allow magically a different tolerance, maybe mesh.tolerance?

Atol = 1e-5;
Rtol = 1e-3;

bool_out = (abs(x-y)<=Atol) || (abs(x-y) <= Rtol*max(abs([x y])));

end