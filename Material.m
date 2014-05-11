classdef Material
    properties
        E           % Elastic Modulus [Pa]
        nu          % Poisson coefficient
        rho         % Density [kg/m3]
        D           % Piezo Coefficients [C/m2] [3x6]
        e           % Electric permitivity [C/Nm2]
    end
    methods 
        function obj = Material(E,nu,rho)
            require(isnumeric([E,nu,rho]), ... 
                'ArgumentError: E,nu, and rho should be numeric')
            require(all([E,nu,rho]>=0), ... 
                'ArgumentError: E,nu, and rho should be greater than 0')
            obj.E = E;
            obj.nu = nu;
            obj.rho = rho;
        end
    end
    methods (Static)
        function obj = Piezo(E,nu,rho,D,permitivity)
            require(isnumeric(D), ...
                'ArgumentError: D should be numeric')
            require(all(size(D)==[3,6]), ...
                'ArgumentError: D should be size [3x6]')
            require(isnumeric(permitivity), ...
                'ArgumentError: permitivity should be numeric')
            obj = Material(E,nu,rho);
            obj.D = D;
            obj.e = permitivity;
        end
    end
end