classdef FemCaseTest < matlab.unittest.TestCase
    methods (Test)
        function BeamTest(testCase)
            %% PRELIMINARY ANALYSIS AND PARAMETERS 
            % Beam Theory analysis
            a = 2;          % X Side length [m]
            b = 0.5;       % Y Side length [m]
            t = 0.1;        % Shell Thickness [m]
            E = 2e7;        % Elasticity Modulus [Pa]
            nu = 0;       % Poisson Coefficient
            rho = 1;        % Density [kg/m3]
            area = b*t;     % Area where Distributed Force is applied
            I = b*(t^3)/12; % Moment of Inertia of the loaded face
            
            % Load is a moment in the end
            tau_z = 1.5;    % Distributed Force [N/m2]
            M = area*tau_z;
            % Expected Displacemente along the x and y axis
            expected_dz     = M*(a^2)/(2*E*I);
            expected_phi    = -M*(a)/(E*I);
            
            %% FEM and MESH
            % Elements along the side
            dofs_per_node = 5;
            dofs_per_ele = 0;
            n = 1;
            m = 2*n;
            mesh = Factory.ShellMesh('AHMAD8',[m,n],[a,b,t]);
            material = Material(E,nu,rho);
            K = @(element) Physics.K_Shell(element,material,3);
            physics = Physics(5,0,K);
            fem = FemCase(mesh,physics);
            
            tol = 1e-5;
            
            %% BC            
            % Fixed End
            x0_edge = (@(x,y,z) (abs(x) < tol));
            base = mesh.find_nodes(x0_edge);
            fem.bc.node_vals.set_val(base,true);
            
            %% LOADS
            % Load the other side
            x1_edge = (@(x,y,z) (abs(x-a) < tol));
            border = find(mesh.find_nodes(x1_edge));
            q = [0 0 0 0 -tau_z]';
            load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
                                                        element,2,q,sc,sv);
            L = mesh.integral_along_surface(dofs_per_node,dofs_per_ele, ...
                                                        border,load_fun);
            fem.loads.node_vals.dof_list_in(L);
            
            %% SOLVE
            fem.solve();    % SIDE EFFECTS
            
            %% TEST
            % Check if the obtained values are the expected ones
            expected_dz
            expected_phi
            fem.dis.node_vals.vals;
            max_dis = max(fem.dis.all_dofs)
            min_dis = min(fem.dis.all_dofs)
            testCase.verifyEqual(true,near(expected_dz,max_dis));
            testCase.verifyEqual(true,near(expected_phi,min_dis));
            
%             %% PLOT
%             PlotMesh(mesh.coords + 1000*fem.dis.node_vals.vals(:,1:3), ...
%                 mesh.connect, ...
%                 ~fem.bc.node_vals.vals, ...
%                             fem.reactions.node_vals.vals);
        end
        function PatchTest(testCase)
            %% PRELIMINARY ANALYSIS AND PARAMETERS 
            % In-plane patch test for one element
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
            sigma = q;

            %% FEM and MESH
            dofs_per_node = 5;
            dofs_per_ele = 0;
            m = 2;
            mesh = Factory.ShellMesh('Q4',[m,m],[a,b,t]);
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
            physics = Physics(dofs_per_node,dofs_per_ele,K);
            fem = FemCase(mesh,physics);
            
            %% BC
            % Simply Supported
            x0_edge = (@(x,y,z) (abs(x) < tol));
            base = mesh.find_nodes(x0_edge);
            fem.bc.node_vals.set_val(base,[true,false,false,false,false]);
            % Fix one node
            x0_point = (@(x,y,z) (norm([x,y]) < tol));
            corner = mesh.find_nodes(x0_point);
            fem.bc.node_vals.set_val(corner,true(5,1));
            %% LOADS
            % Load the other side            
            x1_edge = (@(x,y,z) (abs(x-a) < tol));
            border = find(mesh.find_nodes(x1_edge));
            q = [sigma 0 0 0 0]';
            load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
                element,2,q,sc,sv);
            L = mesh.integral_along_surface(dofs_per_node,dofs_per_ele, ...
                border,load_fun);
            fem.loads.node_vals.dof_list_in(L);
            %% SOLUTION
            fem.solve();    % SIDE EFFECTS
            %% TEST
            % Check if the obtained values are the expected ones
            max_dis = max(fem.dis.all_dofs);
            min_dis = min(fem.dis.all_dofs);
            testCase.verifyEqual(true,near(expected_dx,max_dis));
            testCase.verifyEqual(true,near(expected_dy,min_dis));
        end
    end
end