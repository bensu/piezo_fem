function w = beam_shapes(L,n,x)
    % w = beam_shapes(n,x)
    % Computes one modal shape - n -  of a clamped beam
    lambda = [1.875104069,4.694091133,7.854757438,10.99554073,14.13716839]';
    x_x = lambda(n)*x/L;
    x_L = lambda(n);
    w = (cos(x_x) - cosh(x_x)) - (cos(x_L) + cosh(x_L))* ...
                (sin(x_x) - sinh(x_x)) / (sin(x_L) + sinh(x_L));                                    
    [w_max,loc] = max(abs(w));
    w = sign(w(loc))* w / w_max;
end