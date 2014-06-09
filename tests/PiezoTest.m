classdef PiezoTest < matlab.unittest.TestCase
    methods (Test)
        function BimorphTest(testCase)
            %% Parameters
            a = 100e-3;             % [m] Beam's length
            b = 5e-3;               % [m] Beam's width
            t = 1e-3;               % [m] Beam's thickness
            E   = 2e9;              % [Pa] Elasticity Coefficient
            nu  = 0;                % Poisson coefficient
            rho = 1;                % [kg/m3] Density
            d13 = -0.046*(1e3);     % [(1e-3)C/m2] Piezo constant
            e3  = (1e3)*0.01062e-9; % [(1e-3)C/Nm2] Permitivity
            V = 1e-3;
            %% Expected Values
            % [Development of a Light-Weight Robot End-Effector using 
            % Polymeric Piezoelectric Bimorph, Tzou, 1989][Eq 20]
            expected_w = -3*d13*V*(a^2)/(2*E*t^2); 
            %% Laminate and Material
            piezo_matrix = zeros(3,6);
            piezo_matrix(3,1:3) = [d13 d13 0];
            pzt = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);
            l = 0.5;
            laminate = Laminate([pzt, pzt],[l*t,(1-l)*t]);
            %% Geometry and Mesh
            mesh = Factory.ShellMesh(EleType.AHMAD8,laminate,[1,1],[a,b,t]);
            mesh.laminate_list = laminate;
            %% Physics and FEM
            dofs_per_node = 5;
            dofs_per_ele  = laminate.n_layers;
            gauss_order   = 2;
            K = @(element) Physics.K_PiezoShell(element,gauss_order);
            physics = Physics(dofs_per_node,dofs_per_ele,K);
            fem = FemCase(mesh,physics);
            %% BC
            tol = 1e-9;
            % Clamped
            f_edge = (@(x,y,z) (abs(x) < tol));
            base = mesh.find_nodes(f_edge);
            fem.bc.node_vals.set_val(base,true);
            %% Loads - Initial Displacements
            % All elements charged, initial condition
            fem.initial_dis.ele_vals.vals(:,1) = V/2;
            fem.initial_dis.ele_vals.vals(:,2) = -V/2;
            %% Solve
            fem.solve
            %% Verify
            expected_w;
            w_max = max(abs(fem.dis.node_vals.vals(:,3)));
            testCase.verifyTrue(near(expected_w,w_max))
        end
        function InPlanePatchTest(testCase)
            % Uniform traction applied at edge x = a
            %% Parameters
            s = 1e8;        % [Pa] Applied and expected tension
            E = 123e9;      % [Pa] Elastic Modulus
            nu = 0.3;       % Poisson coefficient
            rho = 1;        % [kg/m3] Density
            t = 0.01;       % [m] Thickness
            a = 0.24;       % [m] Plate length x-side
            b = 0.1;        % [m] Plate length y-side
            d13 = -5*(1e3);         % Piezo constant [(1e-3)C/m2]
            e3 = (12.5e-9)*(1e3);   % Permitivity constant [(1e-3)C/Nm2]
            
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
            piezo_matrix = zeros(3,6);
            piezo_matrix(3,1:3) = [d13 d13 0];
            material = Material.Piezo(E,nu,rho,piezo_matrix,[0 0 e3]);
            laminate = Laminate(material,t);
            mesh = Factory.ShellMesh(EleType.AHMAD4,laminate,[4,2],[a,b,t]);
            % Create the FemCase
            dofs_per_node = 5;
            dofs_per_ele = 1;
            K = @(element) Physics.K_PiezoShell(element,2);
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
            
            %% Check Stresses and Strains (for Report)
            ele_id = 1;
            ele = mesh.ele(ele_id);
            dofs = mesh.all_eles_dofs(dofs_per_node,dofs_per_ele,ele_id);
            D = fem.dis.all_dofs;
            ele_dis = D(dofs);
            B = Physics.B_PiezoShell(ele,laminate.n_layers,1,0,0,0);
            C = laminate.PiezoMatrix(0);
            format long
            x;
            B*ele_dis;
            
            %% Compare
            testCase.verifyTrue(near(expected_u,u_max))
            testCase.verifyTrue(near(expected_v,v_min))
            testCase.verifyTrue(near(expected_V,V))
        end
    end
end