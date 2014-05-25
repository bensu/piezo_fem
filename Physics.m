classdef Physics
    properties
        dofs_per_node
        dofs_per_ele
        k
        m
    end
    methods
        function obj = Physics(dofs_per_node,dofs_per_ele,fun_in)
            require(all(~mod([dofs_per_node,dofs_per_ele],1)), ...
                'ArgumentError: dof numbers should be integers');
            require(isa(fun_in,'function_handle'), ...
                'ArgumentError: fun_in should be a function handle');
            obj.dofs_per_node = dofs_per_node;
            obj.dofs_per_ele = dofs_per_ele;
            obj.k = fun_in;
        end
    end
    methods (Static)
        function obj = Dynamic(dofs_per_node,dofs_per_ele,k,m)
            require(isa(m,'function_handle'), ...
                'ArgumentError: m should be a function handle');
            obj = Physics(dofs_per_node,dofs_per_ele,k);
            obj.m = m;            
        end
        function K = K_PiezoShell(element,material,order)
            % K = K_PiezoShell(element,material,order)
            % K [ele_dof x ele_dof][Float] Stiffness as calculated in
            % Cook 361 12.4-14
            % element [Element]: Requires methods jacobian and B
            % material[Material]: Requires methods E, nu, and D
            % order [Int]: Gauss integration order
            % Constitutive Relationship
            % Both material properties skip 3rd col because it is a plain
            % stress problem
            elastic = Physics.ElasticShell(material);
            piezo = -material.D(:,[1 2 4 5 6]);  
            electric = -diag(material.e);
            C = [   elastic  piezo';
                    piezo   electric];
            % Function to be integrated
            function K_in_point = K_in_point(ksi,eta,zeta)
                % Piezo Part, generates the Electric Field (only z
                % component from 2 or 1 voltage element dof.
                jac = element.jacobian(ksi,eta,zeta);
                cosines = Element.direction_cosines(jac);
                inv_jac = jac \ eye(3);
                dof_per_ele = 1;
                if dof_per_ele == 2
                    dN_ele = zeros(3,2);
                    % V_bottom is first, V_top goes second.
                    % Together they form E_z = V_top - V_bottom
                    dN_ele(3,:) = [-1 1];
                else
                    dN_ele = [0 0 1]';
                    dN_xyz = inv_jac*dN_ele;
                end
                % Mechanics Part
                B_mech = Physics.B_Shell(element,ksi,eta,zeta); % Cook [7.3-10]
                % Join both and trasform the coordinates
                B = blkdiag(Element.T(cosines)*B_mech,cosines*dN_xyz);
                K_in_point = B'*C*B*det(jac);
            end
            fun_in = @(xi,eta,mu) (K_in_point(xi,eta,mu));
            K = Integral.Volume3D(fun_in,order,[-1 1]);
        end
        function L = apply_surface_load(element,order,q,s_coord,s_val)
            % L = apply_load(element,order,q)
            % Generates a load dof vector by integrating a constant load along
            % an element's surface.
            % L [n_ele_dofs x 1][Float]: Load vector
            % element [Element]
            % order [Int]: Gauss integration order
            % q [dof x 1][Float]: Constant applied load
            % s_coord [Int]: coordinate that defines the surface, if ksi,2
            % s_val [Float]: value that the coordinate takes, i.e. ksi = -1;
            function L = apply_point_load(element,q,s_coord,ksi,eta,zeta)
                % L = apply_point_load(element,q,ksi,eta,zeta)
                % Used as lambda in apply_surface_load and apply_volume_load
                jac = element.jacobian(ksi,eta,zeta);
                aux = 1:3;  aux(s_coord) = [];
                v1 = jac(aux(1),:)';
                v2 = jac(aux(2),:)';
                NN = Element.shape_to_diag(length(q),element.N(ksi,eta));
                L = NN'*norm(cross(v1,v2))*q;
            end
            switch (s_coord)
                case 1
                    ksi = s_val;
                    fun_in = @(eta,zeta) (apply_point_load( ...
                                        element,q,s_coord,ksi,eta,zeta));
                case 2
                    eta = s_val;
                    fun_in = @(ksi,zeta) (apply_point_load( ...
                                        element,q,s_coord,ksi,eta,zeta));
                case 3
                    zeta = s_val;
                    fun_in = @(ksi,eta) (apply_point_load( ...
                                        element,q,s_coord,ksi,eta,zeta));
            end
            L = Integral.Surface2D(fun_in,order,[-1 1]);
        end
        function L = apply_volume_load(element,order,q)
            % L = apply_load(element,order,q)
            % Generates a load dof vector by integrating a constant load along
            % the element.
            % L [n_ele_dofs x 1][Float]: Load vector
            % element [Element]
            % order [Int]: Gauss integration order
            % q [dof x 1][Float]: Constant applied load
            function L = apply_point_load(element,q,ksi,eta,zeta)
                % L = apply_point_load(element,q,ksi,eta,zeta)
                % Used as lambda in apply_surface_load and apply_volume_load
                NN = Element.shape_to_diag(length(q),element.N(ksi,eta));
                L = NN'*det(element.jacobian(ksi,eta,zeta))*q;
            end
            fun_in = @(ksi,eta,zeta) (apply_point_load(element,q, ...
                ksi,eta,zeta));
            L = Integral.Volume3D(fun_in,order,[-1 1]);
        end
        function M = M_Shell(element,material,order)
            % M = M_Shell(element,material,order)
            % M [n_dofxn_dof][Float] Mass as calculated in Cook 361 13.2-5
            % element [Element]: Requires methods jacobian and B
            % material[Material]: Requires property rho
            % order [Int]: Gauss integration order
            rho = material.rho;
            function M_in_point = M_in_point(ksi,eta,zeta)
                jac = element.jacobian(ksi,eta,zeta);
                N = element.ShellN(ksi,eta,zeta);
                M_in_point = rho*(N')*N*det(jac);
            end
            fun_in = @(xi,eta,mu) (M_in_point(xi,eta,mu));
            M = Integral.Volume3D(fun_in,order,[-1 1]);
        end
        function K = K_Shell(element,material,order)
            % K = K_Shell(element,material,order)
            % K [n_dofxn_dof][Float] Stiffness as calculated in Cook 361 12.4-14
            % element [Element]: Requires methods jacobian and B
            % material[Material]: Requires properties E and nu
            % order [Int]: Gauss integration order
            C = Physics.ElasticShell(material);
            function K_in_point = K_in_point(ksi,eta,zeta)
                jac = element.jacobian(ksi,eta,zeta);
                cosines = Element.direction_cosines(jac);
                B = Element.T(cosines)* ...
                    Physics.B_Shell(element,ksi,eta,zeta); % Cook [7.3-10]
                K_in_point = B'*C*B*det(jac);
            end
            fun_in = @(xi,eta,mu) (K_in_point(xi,eta,mu));
            K = Integral.Volume3D(fun_in,order,[-1 1]);
        end
        function K = K_Shell_selective(element,material,normal_order,shear_order)
            % K = K_Shell_selective(element,material,normal_order,shear_order)
            % K [n_dofxn_dof][Float] Stiffness as calculated in Cook 361 12.4-14
            % with selective integration
            % element [Element]: Requires methods jacobian and B
            % material[Material]: Requires methods E and nu
            % normal_order [Int]: Gauss integration order for normal part
            % shear_order  [Int]: Gauss integration order for shear part
            C = Physics.ElasticShell(material);
            function K_in_point = K_in_point(c_matrix,ksi,eta,zeta)
                jac = element.jacobian(ksi,eta,zeta);
                cosines = Element.direction_cosines(jac);
                B = Element.T(cosines)* ...
                    Physics.B_Shell(element,ksi,eta,zeta); % Cook [7.3-10]
                K_in_point = B'*c_matrix*B*det(jac);
            end
            % Shear part
            C_shear = C;	C_shear(1:3,1:3) = 0;
            fun_shear = @(xi,eta,mu) (K_in_point(C_shear,xi,eta,mu));
            K_s = Integral.Volume3D(fun_shear,shear_order,[-1 1]);
            % Transverse part
            C_n = C;    C_n(4:5,4:5) = 0;
            fun_n = @(xi,eta,mu) (K_in_point(C_n,xi,eta,mu));
            K_n = Integral.Volume3D(fun_n,normal_order,[-1 1]);
            K = K_s + K_n;
        end
        function C = Elastic(material)
            % Computes the Elastic Tensor in matrix form for an Isotropic
            % material
            E = material.E;
            nu = material.nu;
            C = E/((1+nu)*(1-2*nu))* ...
                    [1-nu   nu      nu      0          0          0
                    nu      1-nu	nu      0          0          0
                    nu      nu      1-nu    0          0          0
                    0       0       0       (1-2*nu)/2 0          0
                    0       0       0       0          (1-2*nu)/2 0
                    0       0       0       0          0          (1-2*nu)/2];
        end
        function C = ElasticPlainStress(material)
            % Computes the Elastic Tensor in matrix form for an Isotropic
            % material
            E = material.E;
            nu = material.nu;
            aux = E/(1-nu^2)* ...
                    [1  nu  nu
                     nu 1   nu
                     nu nu  1];
            C = blkdiag(aux,0.5*E*(1+nu)*eye(3));
        end
        function C = ElasticShell(material)
            % Returns the Elastic Tensor for Plain Stress in a shell
            % Cook pg 361 12.5-12
            C = Physics.ElasticPlainStress(material);
            c = 5/6;
            C(5,5) = c*C(5,5);
            C(6,6) = c*C(6,6);
            C(3,:) = [];
            C(:,3) = [];
        end
        function B = B_Shell(element,ksi,eta,zeta)
            % B = B_Shell(element,ksi,eta,zeta)
            % B [Float][6 x n_ele_dofs]: Relates the element's dof values 
            % with the mechanical strain vector. Cook [12.5-10]
            % element [Element]
            % ksi, eta, zeta [Float][Scalar] in [-1,1] local coordinates            
            dofs_per_node = 5;
            % Prepare values
            v = element.normals;
            N  = element.N(ksi,eta);
            invJac = element.jacobian(ksi,eta,zeta) \ eye(3);
            dN = invJac(:,1:2)*element.dN(ksi,eta);
            B  = zeros(6,element.n_nodes*dofs_per_node);
            % B matrix has the same structure for each node, 
            % written as [aux1 aux2].
            % Loop through the mesh.connect coords and get each B_node, 
            % then add it to its columns in the B matrix
            for n = 1:element.n_nodes
                v1 = v(:,1,n);  % In Cook [12.5-3] as {l1i,m1i,n1i} 
                v2 = v(:,2,n);  % In Cook [12.5-3] as {l2i,m2i,n2i} 
                dZN = dN(:,n)*zeta + N(n)*invJac(:,3);
                % aux1: Part of node's B unrelated to rotational dofs and zeta
                aux1 = [ dN(1,n)    0           0
                         0          dN(2,n)     0
                         0          0           dN(3,n)
                         dN(2,n)    dN(1,n)     0
                         0          dN(3,n)     dN(2,n)
                         dN(3,n)    0           dN(1,n) ];
                % aux2: Part of node's B related to rotational dofs and zeta
                aux2 = [ -v2.*dZN                        v1.*dZN
                         -v2(1)*dZN(2) - v2(2)*dZN(1)    v1(1)*dZN(2) + v1(2)*dZN(1)
                         -v2(2)*dZN(3) - v2(3)*dZN(2)    v1(2)*dZN(3) + v1(3)*dZN(2)
                         -v2(1)*dZN(3) - v2(3)*dZN(1)    v1(1)*dZN(3) + v1(3)*dZN(1) ]*0.5*element.thickness(n);
                % Add that node's part to the complete B
                B(:,index_range(dofs_per_node,n)) = [aux1 aux2];
            end
        end
        function H_out = H
            % H_out = H
            % H_out [6x9] Cook pg 181 [6.7-5]
            % goes from diff_U_xyz to Strain vector
            H_out = zeros(6,9);
            H_out([1 10 18 22 26 35 42 47 51]) = 1;
        end
    end
end
    