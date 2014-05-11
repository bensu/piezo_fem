classdef FemCaseTest < matlab.unittest.TestCase
    methods (Test)
        function BeamTest(testCase)
            % Patch test for one element
            a = 2;          % X Side length [m]
            b = 0.4;        % Y Side length [m]
            t = 0.1;        % Shell Thickness [m]
            q = 1.5;        % Distributed Force [N/m2]
            E = 2e7;          % Elasticity Modulus [Pa]
            nu = 0.3;       % Poisson Coefficient
            rho = 1;        % Density [kg/m3]
            area = b*t;     % Area where Distributed Force is applied
            I = b*(t^3)/12; % Moment of Inertia of the loaded face
            % Expected Displacemente along the x and y axis
            % http://en.wikipedia.org/wiki/Deflection_%28engineering%29
            %#End_loaded_cantilever_beams
            expected_dz = (q*area)*(a^3)/(3*E*I)
            expected_phi = -(q*area)*(a^2)/(2*E*I)
            
            % Elements along the side
            m = 100;
            n = 40;
            % With 
            F = (q*area)/(n*2);   % Force at each node.

            % Create the FemCase
            mesh = Factory.ShellMesh([m,n],[a,b,t]);
            material = Material(E,nu,rho);
            K = @(element) Physics.K_Shell(element,material,2);
            physics = Physics(5,0,K);
            fem = FemCase(mesh,physics);
            
            tol = 1e-5;
            % Fixed End
            x0_edge = (@(x,y,z) (abs(x) < tol));
            base = mesh.find_nodes(x0_edge);
            fem.bc.node_vals.set_val(base,true);
            % Load the other side
            x1_edge = (@(x,y,z) (abs(x-a) < tol));
            border = mesh.find_nodes(x1_edge);
            fem.loads.node_vals.set_val(border,[0,0,F,0,0]);
            % Load consistently
            x1_inner_edge = @(x,y,z) ((abs(x-a) < tol) && ...
                                     ~(abs(y) < tol)   && ...
                                     ~(abs(y-b) < tol));
            inner_border = mesh.find_nodes(x1_inner_edge);
            fem.loads.node_vals.set_val(inner_border,[0,0,2*F,0,0]);
            % Solve
            fem.solve();    % SIDE EFFECTS
            % Check if the obtained values are the expected ones
            PlotMesh(mesh.coords + 1000*fem.dis.node_vals.vals(:,1:3), ...
                mesh.connect, ...
                ~fem.bc.node_vals.vals, ...
                            fem.reactions.node_vals.vals);
            fem.dis.node_vals.vals;
            max_dis = max(fem.dis.all_dofs)
            min_dis = min(fem.dis.all_dofs)
            testCase.verifyEqual(true,near(expected_dz,max_dis));
            testCase.verifyEqual(true,near(expected_phi,min_dis));
        end
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
            expected_dx = q*a/E;
            expected_dy = -nu*b*q/E;
            
            % Elements along the side
            m = 2;
            % With 
            F = q*area/(m*2);   % Force at each node.

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
            fem.bc.node_vals.set_val(base,[true,false,false,false,false]);
            % Fix one node
            x0_point = (@(x,y,z) (norm([x,y]) < tol));
            corner = mesh.find_nodes(x0_point);
            fem.bc.node_vals.set_val(corner,true(5,1));
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
            % Solve
            fem.solve();    % SIDE EFFECTS
            % Check if the obtained values are the expected ones
            max_dis = max(fem.dis.all_dofs);
            min_dis = min(fem.dis.all_dofs);
            testCase.verifyEqual(true,near(expected_dx,max_dis));
            testCase.verifyEqual(true,near(expected_dy,min_dis));
        end
    end
end