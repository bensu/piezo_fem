classdef LaminateTest < matlab.unittest.TestCase
    methods (Test)
        function BimorphTest(testCase)
            %% Geometry and Mesh
            a = 100e-3;             % [m] Beam's length
            b = 5e-3;               % [m] Beam's width
            t = 1e-3;               % [m] Beam's thickness
            mesh = Factory.ShellMesh(EleType.AHMAD8,[4,2],[a,b,t]);
            
            %% Laminate and Material
            E = 2e9;                % [Pa] Elasticity Coefficient
            nu = 0;                 % Poisson coefficient
            d13 = -0.046*(1e3);     % [(1e-3)C/m^2] Piezo coefficient
            e3 = (1e3)*(0.1062e-9); % [(1e-3)C/Nm^2] Permitivity constant
            piezo_matrix = zeros(3,6);
            piezo_matrix(3,1:3) = [d13 d13 0];
            piezo = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);
            laminate = Laminate([piezo, piezo],[t/2,t/2]); % is it /2 ???
            
            %% Physics and FEM
            dofs_per_node = 5;
            dofs_per_ele = 2;
            K = @(element) Physics.K_PiezoShell(element,laminate,2);
            physics = Physics(dofs_per_node,dofs_per_ele,K);
            fem = FemCase(mesh,physics);
        end
    end
end