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
            expected_dz = 0.8*(a^2); % [m]
            expected_dy = 0;         % [m] because poisson = 0
            expected_V  = 1.686e3;   % [V]
            
            dofs_per_node = 5;
            dofs_per_ele  = 1;
            
            %% MESH
           
            mesh = Factory.ShellMesh('AHMAD4',[2,1],[a,b,t]);
            
            %% MATERIAL
            piezo_matrix = zeros(3,6);
            piezo_matrix(3,1:3) = [e13 e13 0];
            material = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);

            %% FEM and MESH
            % Create the FemCase
            K = @(element) Physics.K_PiezoShell(element,material,2);
            physics = Physics(dofs_per_node,dofs_per_ele,K);
            fem = FemCase(mesh,physics);
            
            %% BC
            % Clamped since there is no poisson effect
            tol = 1e-9;
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
            
            dz = max(fem.dis.node_vals.vals(:,3));
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
        % Uniform tension applied at edge x = a
            
            %% Parameters
            s = 1e8;        % [Pa] Applied and expected tension
            E = 123e9;      % [Pa] Elastic Modulus
            nu = 0.3;       % Poisson coefficient
            rho = 1;        % [kg/m3] Density
            t = 0.01;       % [m] Thickness
            a = 0.24;       % [m] Plate length x-side
            b = 0.1;        % [m] Plate length y-side
            d13 = -5;       % [???] Piezo Constant
            e3 = 12.5e-9;    % [???] Electric permitivity
            
            %% Expected Values
            % Solve constitutive equation by assuming sigma_yy = 0; D_z = 0
            % Constitutive Matrix
            aux = E/(1-nu^2)*[1   nu
                              nu  1];
            C = blkdiag(aux,-e3);
            C(3,1:2) = -d13;
            C(1:2,3) = -d13;
            S = [s 0 0]';   % Applied tension and D_z
            x = C \ S;  % Solution
            % Solution decomposition
            expected_u = x(1)*a;
            expected_v = x(2)*b;
            expected_V = x(3)*t/2;  % Not clear why is the /2 needed
            
            %% FEM and MESH
            mesh = Factory.ShellMesh('AHMAD4',[4,2],[a,b,t]);
            piezo_matrix = zeros(3,6);
            piezo_matrix(3,1:3) = [d13 d13 0];
            material = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);
            % Create the FemCase
            dofs_per_node = 5;
            dofs_per_ele = 1;
            K = @(element) Physics.K_PiezoShell(element,material,2);
            physics = Physics(dofs_per_node,dofs_per_ele,K);
            fem = FemCase(mesh,physics);
            
            %% BC
            % Simply supported
            tol = 1e-9;
            % One edge restricted on x direction
            f_edge = (@(x,y,z) (abs(x) < tol));
            base = mesh.find_nodes(f_edge);
            fem.bc.node_vals.set_val(base,[true false false false false]);
            % One fixed point
            f_corner = (@(x,y,z) (norm([x y z]) < tol));
            corner = mesh.find_nodes(f_corner);
            fem.bc.node_vals.set_val(corner,true);
            
            %% LOADS
            % Load the other side
            x1_edge = (@(x,y,z) (abs(x-a) < tol));
            border = find(mesh.find_nodes(x1_edge));
            s = [s 0 0 0 0]';
            load_fun = @(element,sc,sv) Physics.apply_surface_load( ...
                element,2,s,sc,sv);
            L = mesh.integral_along_surface(dofs_per_node,0, ...
                border,load_fun);
            fem.loads.node_vals.dof_list_in(L);
            
            %% Solve   
            fem.solve();
            u_max = max(fem.dis.node_vals.vals(:,1));
            v_min = min(fem.dis.node_vals.vals(:,2));
            V = fem.dis.ele_vals.vals(1,1);
            
            %% Compare
            testCase.verifyTrue(near(expected_u,u_max))
            testCase.verifyTrue(near(expected_v,v_min))
            testCase.verifyTrue(near(expected_V,V))
        end
    end
end