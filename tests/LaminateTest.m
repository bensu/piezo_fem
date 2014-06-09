classdef LaminateTest < matlab.unittest.TestCase
    methods (Test)
        function BasicTest(testCase)
            %% Laminate and Material
            a = 100e-3;             % [m] Beam's length
            b = 5e-3;               % [m] Beam's width
            t = 1e-3;               % [m] Beam's thickness
            E = 2e9;                % [Pa] Elasticity Coefficient
            nu = 0;                 % Poisson coefficient
            rho = 1;                % [kg/m3] Density
            metal = Material(E,nu,rho);
            laminate = Laminate([metal, metal],[t/2,t/2]); % is it /2 ???
            %% Mesh
            mesh = Factory.ShellMesh(EleType.AHMAD8,laminate,[4,2],[a,b,t]);
            %% Physics and FEM
            dofs_per_node = 5;
            dofs_per_ele = 0;
            K = @(element) Physics.K_Shell(element,2);
            physics = Physics(dofs_per_node,dofs_per_ele,K);
            fem = FemCase(mesh,physics);
        end
    end
end