classdef FemCaseTest < matlab.unittest.TestCase
    methods (Test)
        function PatchTest(testCase)
            % Patch test for one element
            a = 0.3;        % X Side length [m]
            b = 0.4;        % Y Side length [m]
            t = 0.1;        % Shell Thickness [m]
            q = 1.5;        % Distributed Force [N/m2]
            E = 2;          % Elasticity Modulus [Pa]
            nu = 0.3;       % Poisson Coefficient
            rho = 1;        % Density [kg/m3]
            area = b*t;     % Area where Distributed Force is applied
            % Expected Displacemente along the x and y axis
            expected_dx = (1-nu^2)*q*a/E;
            expected_dy = -nu*expected_dx/(1-nu)*b/a;
            
            % Elements along the side
            m = 2;
            % With 
            F = q*area/(m*2);   % Force at each element.

            mesh = Factory.ShellMesh([m,m],[a,b,t]);
            % Break the symmetry
            tol = 1e-5;
            inner = @(x,y,z) (~(abs(x-a) < tol) && ~(abs(x) < tol) && ...
                              ~(abs(y-b) < tol) && ~(abs(y) < tol));
            inner_nodes = mesh.find_nodes(inner);
            % THIS MIGHT BREAK IF THE RESULTING NODE COORD ENDS UP IN A
            % STRANGE PLACE BECAUSE OF THE 1.1
            mesh.coords(inner_nodes(1),:) = mesh.coords(inner_nodes(1),:)*1.1;
            material = Material(E,nu,rho);
            
            % Create the FemCase
            K = @(element) Physics.K_Shell(element,material,2);
            physics = Physics(5,0,K);
            fem = FemCase(mesh,physics);
            
            % Simply Supported
            x0_edge = (@(x,y,z) (abs(x) < tol));
            base = mesh.find_nodes(x0_edge);
            fem.bc.node_vals.set_val(base,[false,true,true,true,true]);
            % Fix one node
            x0_point = (@(x,y,z) (norm([x,y]) < tol));
            corner = mesh.find_nodes(x0_point);
            fem.bc.node_vals.set_val(corner,false(5,1));
            % Load the other side
            x1_edge = (@(x,y,z) (abs(x-a) < tol));
            border = mesh.find_nodes(x1_edge);
            fem.loads.node_vals.set_val(border,[F,0,0,0,0]);
            % Load consistently
            x1_inner_edge = @(x,y,z) ((abs(x-a) < tol) && ...
                                     ~(abs(y) < tol)   && ...
                                     ~(abs(y-b) < tol));
            inner_border = mesh.find_nodes(x1_inner_edge);
            fem.loads.node_vals.set_val(inner_border,[2*F,0,0,0,0]);
            % BC and Loads to DOF form
            L = fem.loads.all_dofs;
            F = fem.bc.all_dofs;
            % Create Stiffness
            S = fem.mesh.assembly(fem.physics.dofs_per_node, ...
                                  fem.physics.dofs_per_ele,  ...
                                  fem.physics.k);
            D = zeros(size(S,1),1);
            D(F,1) = S(F,F) \ L(F);
            fem.dis.dof_list_in(D);
            % Check if the obtained values are the expected ones
            max_dis = max(fem.dis.all_dofs);
            min_dis = min(fem.dis.all_dofs);
            testCase.verifyEqual(true,near(expected_dx,max_dis));
            testCase.verifyEqual(true,near(expected_dy,min_dis));
        end
    end
end