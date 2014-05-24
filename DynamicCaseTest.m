classdef DynamicCaseTest < matlab.unittest.TestCase
    methods (Test)
        function MassMatrix(testCase)
            %% PRELIMINARY ANALYSIS AND PARAMETERS 
            % Beam Theory analysis
            a = 2;          % X Side length [m]
            b = 0.5;       % Y Side length [m]
            t = 0.1;        % Shell Thickness [m]
            E = 2e7;        % Elasticity Modulus [Pa]
            nu = 0;       % Poisson Coefficient
            rho = 1;        % Density [kg/m3]
         
            %% FEM and MESH
            % Elements along the side
            dofs_per_node = 5;
            dofs_per_ele = 0;
            n = 1;
            m = 2*n;
            mesh = Factory.ShellMesh('AHMAD8',[m,n],[a,b,t]);
            material = Material(E,nu,rho);
            M = @(element) Physics.M_Shell(element,material,3);
            K = @(element) Physics.K_Shell(element,material,3);
            physics = Physics.Dynamic(dofs_per_node,dofs_per_ele,K,M);
            fem = FemCase(mesh,physics);
            
        end
    end
end