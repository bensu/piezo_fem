classdef Integral
    methods (Static)
        function [w, gp, n] = gauss(ng)
            % [w, gp, n] = gauss(ng)
            % ng: [Int][dim x 1] where each component indicates how many gauss points
            % shold be used to integrate along that direction
            % wgauss [n x 1]:     Gauss weights
            % gpts [n x dim]:     Evaluation points for each coordinates
            % n [Int]:            Number of evaluation points
            % The function works by calling Integral.gauss1D for each dimension and then
            % assembling w and gp.
            % Note: n is redundant as an output since it is contained in length(wgauss)
            dims = length(ng);
            n = prod(ng);
            w = ones(n,1);
            gp = zeros(n,dims);
            switch dims
                case 3
                    [w1,gp1] = Integral.gauss1D(ng(1));
                    [w2,gp2] = Integral.gauss1D(ng(2));
                    [w3,gp3] = Integral.gauss1D(ng(3));
                    
                    counter = 1;
                    for ig1 = 1:ng(1)
                        for ig2 = 1:ng(2)
                            for ig3 = 1:ng(3)
                                w(counter) = w(counter)*w1(ig1)*w2(ig2)*w3(ig3);
                                gp(counter,:) = [gp1(ig1) gp2(ig2) gp3(ig3)];
                                counter = counter + 1;
                            end
                        end
                    end
                case 2
                    [w1,gp1] = Integral.gauss1D(ng(1));
                    [w2,gp2] = Integral.gauss1D(ng(2));
                    
                    counter = 1;
                    for ig1 = 1:ng(1)
                        for ig2 = 1:ng(2)
                            w(counter) = w(counter)*w1(ig1)*w2(ig2);
                            gp(counter,:) = [gp1(ig1) gp2(ig2)];
                            counter = counter + 1;
                        end
                    end
                case 1
                    [w,gp] = Integral.gauss1D(n);
                    w = w';
                    gp = gp';
            end
        end
        function [w,gp] = gauss1D(n)
            % [w,gp] = Integral.gauss1D(n)
            % n [Int]: number of points along the line
            % w [Float][n x 1]:     Gauss weights for integration along the line
            % gp [Float][n x 1]:    Places where to evaluate the function
            
            % When approximating integral(f(x),x_0,x_1) with sum(w*f(gp(i))) this
            % function returns the values where to evaluate the function (gp) and the
            % weight that should be assigned to each result (w)
            switch n 
                case 1
                    w  = 2;
                    gp = 0;
                case 2
                    w  = [1 1];
                    a  = sqrt(3)/3;
                    gp = [-a a];
                case 3
                    w  = [5/9 8/9 5/9];
                    a  = sqrt(3/5);
                    gp = [-a 0 a];
                case 4
                    a  = sqrt((3 - 2*sqrt(6/5))/7);
                    b  = sqrt((3 + 2*sqrt(6/5))/7);
                    gp = [-b -a a b];
                    wa = (18 + sqrt(30))/36;
                    wb = (18 - sqrt(30))/36;
                    w  = [wb wa wa wb];
                case 5
                    a  = 1/3*sqrt(5 - 2*sqrt(10/7));
                    b  = 1/3*sqrt(5 + 2*sqrt(10/7));
                    gp = [-b -a 0 a b];
                    wa = (322 + 13*sqrt(70))/900;
                    wb = (322 - 13*sqrt(70))/900;
                    w  = [wb wa 128/225 wa wb];
                case 6
                    a  = 0.932469514203152;
                    b  = 0.661209386466265;
                    c  = 0.238619186083197;
                    wa = 0.171324492379170;
                    wb = 0.360761573048139;
                    wc = 0.467913934572691;
                    gp = [-a -b -c c b a];
                    w  = [wa wb wc wc wb wa];
                    %     case 7
                    %         a  = 0.949107912342759;
                    %         b  = 0.741531185599394;
                    %         c  = 0.405845151377397;
                    %         d  = 0.0;
                    %         wa = 0.129484966168870;
                    %         wb = 0.279705391489277;
                    %         wc = 0.381830050505119;
                    %         wd = 0.417959183673469;
                    %         gp = [-a -b -c d c b a];
                    %         w  = [wa wb wc wd wc wb wa];
                    %     case 8
                    %         a  = 0.932469514203152;
                    %         b  = 0.661209386466265;
                    %         c  = 0.238619186083197;
                    %         wa = 0.171324492379170;
                    %         wb = 0.360761573048139;
                    %         wc = 0.467913934572691;
                    %         gp = [-a -b -c c b a];
                    %         w  = [wa wb wc wc wb wa];
                    %     case 9
                    %         a  = 0.932469514203152;
                    %         b  = 0.661209386466265;
                    %         c  = 0.238619186083197;
                    %         wa = 0.171324492379170;
                    %         wb = 0.360761573048139;
                    %         wc = 0.467913934572691;
                    %         gp = [-a -b -c c b a];
                    %         w  = [wa wb wc wc wb wa];
                    %     case 10
                    %         a  = 0.932469514203152;
                    %         b  = 0.661209386466265;
                    %         c  = 0.238619186083197;
                    %         wa = 0.171324492379170;
                    %         wb = 0.360761573048139;
                    %         wc = 0.467913934572691;
                    %         gp = [-a -b -c c b a];
                    %         w  = [wa wb wc wc wb wa];
            end
        end
        function mat_out = Volume3D(fun_in,fun_size,order,interval)
            % mat_out = Volume3D(fun_in,fun_size,order,interval)
            % mat_out [fun_size,fun_size]
            % fun_in [Fhandle] function(xi,eta,mu) to be integrated
            % order [Int][>0] Gauss integration
            % interval [1x2][Float] Cube sides for integration
            require(isnumeric(order), ...
                'ArgumentError: order should be numeric');
            require(order > 0, ...
                'ArgumentError: order should be > 0');
            require(~mod(order,1), ...
                'ArgumentError: order should be integer');
            require(length(interval)==2, ...
                'ArgumentError: interval should have two numbers');
            require(isnumeric(interval), ...
                'ArgumentError: interval should be numeric');
            require(isnumeric(fun_size), ...
                'ArgumentError: fun_size should be numeric');
            require(fun_size > 0, ...
                'ArgumentError: fun_size should be > 0');
            require(~mod(fun_size,1), ...
                'ArgumentError: fun_size should be integer');
            [gauss_p,gauss_w] = Integral.lgwt(order,interval(1),interval(2));     % sampling points & weights
            mat_out = zeros(fun_size); % Initialization of Matrix
            % Numerical integration
            for int_mu = 1:order
                % sampling point in z-axis
                mu = gauss_p(int_mu,1);
                % weight in z-axis
                wtz = gauss_w(int_mu,1);
                for int_xi = 1:order
                    % sampling point in x-axis
                    xi = gauss_p(int_xi,1);
                    % weight in x-axis
                    wtx = gauss_w(int_xi,1);	
                    for int_eta= 1:order
                        % sampling point in y-axis
                        eta = gauss_p(int_eta,1);
                        % weight in y-axis
                        wty = gauss_w(int_eta,1);
                        aux = fun_in(xi,eta,mu);
                        mat_out = mat_out + wtx*wty*wtz*aux;
                    end
                end
            end
            require(all(size(mat_out)==fun_size), ...
                            'ArgumentError: fun_size and fun_in should match')
        end
        function [x,w] = lgwt(N,a,b)
            % lgwt.m
            %
            % This script is for computing definite integrals using Legendre-Gauss
            % Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
            % [a,b] with truncation order N
            %
            % Suppose you have a continuous function f(x) which is defined on [a,b]
            % which you can evaluate at any x in [a,b]. Simply evaluate it at all of
            % the values contained in the x vector to obtain a vector f. Then compute
            % the definite integral using sum(f.*w);
            %
            % Written by Greg von Winckel - 02/25/2004
            N=N-1;
            N1=N+1; N2=N+2;
            xu=linspace(-1,1,N1)';
            % Initial guess
            y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
            % Legendre-Gauss Vandermonde Matrix
            L=zeros(N1,N2);
            % Derivative of LGVM
            Lp=zeros(N1,N2);
            % Compute the zeros of the N+1 Legendre Polynomial
            % using the recursion relation and the Newton-Raphson method
            y0=2;
            % Iterate until new points are uniformly within epsilon of old points
            while max(abs(y-y0))>eps 
                L(:,1)=1;
                Lp(:,1)=0;
                L(:,2)=y;
                Lp(:,2)=1;
                for k=2:N1
                    L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
                end
                Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
                y0=y;
                y=y0-L(:,N2)./Lp;
                
            end
            % Linear map from[-1,1] to [a,b]
            x=(a*(1-y)+b*(1+y))/2;
            % Compute the weights
            w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
        end
    end
end