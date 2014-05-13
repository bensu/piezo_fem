classdef PiezoTest < matlab.unittest.TestCase
    methods (Test)
        function OutofPlanePatchTest(testCase)
            %% PRELIMINARY ANALYSIS AND PARAMETERS
            % In-plane patch test for one element
            % Exactly like 5.1.8 Aachen Thesis
            a = 0.24;       % X Side length [m]
            b = 0.12;       % Y Side length [m]
            t = 0.01;       % Shell Thickness [m]
            M = 1e3;        % Moment applied to points [N]
            area = b*t;
            M_dist = 2*M/(area);
            E	= 123e9;    % Elasticity Modulus [N/m2]
            nu	= 0;        % Poisson Coefficient
            rho = 1;        % Density [kg/m3]
            e13 = -5;       % Piezo constant [C/m2]
            e23 = 0;        % Piezo constant [C/m2], not used for this problem
            e3  = 12.5e-9;  % Permitivity constant [C/Nm2]
            
            % Expected Displacemente along the x and y axis
            expected_dz = 0.8*(a^2)    % [m]
            expected_dy = 0;         % [m] because poisson = 0
            expected_V  = 1.686e3;   % [V]
            
            dofs_per_node = 5;
            dofs_per_ele  = 1;
            
            %% MESH
           
            mesh = Factory.ShellMesh('AHMAD4',[4,2],[a,b,t]);
            
            %% MATERIAL
            piezo_matrix = zeros(3,6);
            piezo_matrix(3,1:3) = [e13 e13 0];
            material = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);

            %% FEM and MESH
            % Create the FemCase
            K = @(element) Physics.K_PiezoShell(element,material,2);
            physics = Physics(dofs_per_node,dofs_per_ele,K);
            fem = FemCase(mesh,physics);
            tol = 1e-9;
            %% BC
            % Clamped since there is no poisson effect
            x0_edge = (@(x,y,z) (abs(x) < tol));
            base = mesh.find_nodes(x0_edge);
            fem.bc.node_vals.set_val(base,true);     
            % All bottom layer is grounded
%             fem.bc.ele_vals.vals(:,1) = true;

            %% LOADS
            % Load the other side
            x1_edge = (@(x,y,z) (abs(x-a) < tol));
            border = find(mesh.find_nodes(x1_edge));
            q = [0 0 0 0 -M_dist]';
            load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
                                                        element,2,q,sc,sv);
            L = mesh.integral_along_surface(dofs_per_node,0, ...
                                                        border,load_fun);
            fem.loads.node_vals.dof_list_in(L);
            
            %% SOLUTION
            fem.solve();    % SIDE EFFECTS
            %% TEST
            % Check if the obtained values are the expected ones
            fem.loads.node_vals.vals
            fem.reactions.node_vals.vals
            
            dz = max(fem.dis.node_vals.vals(:,3))
            V   = (fem.dis.ele_vals.all_dofs);
            %             testCase.verifyEqual(true,near(expected_dz,dz));
%             testCase.verifyEqual(true,near(expected_V,max_V));
            
            %% PLOT
            PlotMesh(mesh.coords + 1*fem.dis.node_vals.vals(:,1:3), ...
                mesh.connect, ...
                ~fem.bc.node_vals.vals, ...
                            fem.reactions.node_vals.vals);
        end
        function InPlanePatchTest(testCase)
            %% PRELIMINARY ANALYSIS AND PARAMETERS
            % In-plane patch test for one element
            % Exactly like 5.1.8 Aachen Thesis
            a = 0.24;       % X Side length [m]
            b = 0.12;       % Y Side length [m]
            t = 0.01;       % Shell Thickness [m]
            F = 6e4;        % Force applied to points [N]
            area = b*t;
            sigma = (2*F)/area;
            E	= 123e9;    % Elasticity Modulus [N/m2]
            nu	= 0;        % Poisson Coefficient
            rho = 1;        % Density [kg/m3]
            e13 = -5;       % Piezo constant [C/m2]
            e23 = 0;        % Piezo constant [C/m2], not used for this problem
            e3  = 12.5e-9;  % Permitivity constant [C/Nm2]
            
            % Expected Displacemente along the x and y axis
            expected_dx = 1.9205e-4;   % [m]
            expected_dy = 0;       % [m] because poisson = 0
            expected_V  = 1.6e3;   % [V]
            
            dofs_per_node = 5;
            dofs_per_ele  = 1;
            
            %% MESH
            % Taken exaclty from example, already asymmetrical
            mesh = Factory.ShellMesh('AHMAD4',[4,2],[a,b,t]);

            %% MESH COPIED FROM BOOK
%             coords = [ 0.04 0.02    0;
%                        0.08 0.08    0;
%                        0.18 0.03    0;
%                        0.16 0.08    0;
%                        0    0       0;
%                        0    b       0;
%                        a    0       0;
%                        a    b       0];
%            
%             connect = [ 5 1 2 6;
%                         5 7 3 1;
%                         1 3 4 2;
%                         6 2 4 8;
%                         7 8 4 3];
%                     
%             n_node = size(coords,1);
%             thickness = t*ones(1,n_node);
%             mesh = Mesh('AHMAD4',coords,connect,thickness);                       
%             
            %% MATERIAL
            piezo_matrix = zeros(3,6);
            piezo_matrix(3,1:3) = [e13 e13 0];
            material = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);

            %% FEM and MESH
            % Create the FemCase
            K = @(element) Physics.K_PiezoShell(element,material,2);
            physics = Physics(dofs_per_node,dofs_per_ele,K);
            fem = FemCase(mesh,physics);
            tol = 1e-9;
            %% BC
            % Clamped since there is no poisson effect
            x0_edge = (@(x,y,z) (abs(x) < tol));
            base = mesh.find_nodes(x0_edge);
            fem.bc.node_vals.set_val(base,true);            
            % All bottom layer is grounded
%             fem.bc.ele_vals.vals(:,1) = true;
            %% LOADS
            % Load the other side
            x1_edge = (@(x,y,z) (abs(x-a) < tol));
            border = find(mesh.find_nodes(x1_edge));
            q = [sigma 0 0 0 0]';
            load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
                                                        element,2,q,sc,sv);
            L = mesh.integral_along_surface(dofs_per_node,0, ...
                                                        border,load_fun);
            fem.loads.node_vals.dof_list_in(L);

            %% SOLUTION
            fem.solve();    % SIDE EFFECTS
            %% TEST
            % Check if the obtained values are the expected ones
            expected_dx;
            max_dis = max(fem.dis.node_vals.all_dofs);
            expected_dy;
            min_dis = min(fem.dis.node_vals.all_dofs);
            expected_V;
            max_V   = max(fem.dis.ele_vals.all_dofs);
            testCase.verifyEqual(true,near(expected_dx,max_dis));
            testCase.verifyEqual(true,near(expected_dy,min_dis));
%             testCase.verifyEqual(true,near(expected_V,max_V));
            
            %% PLOT
%             PlotMesh(mesh.coords + 1000*fem.dis.node_vals.vals(:,1:3), ...
%                 mesh.connect, ...
%                 ~fem.bc.node_vals.vals, ...
%                             fem.reactions.node_vals.vals);
        end
    end
end