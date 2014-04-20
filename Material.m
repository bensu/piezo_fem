classdef Material
    properties
        E           % Elastic Modulus [Pa]
        nu          % Poisson coefficient
        rho         % Density [kg/m3]
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
end