classdef Physics
    properties
        dofs_per_node
        dofs_per_ele
        k
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
        function K = K_PiezoShell(element,material,order)
        % K [ele_dof x ele_dof][Float] Stiffness as calculated in 
        % Cook 361 12.4-14
        % element [Element]: Requires methods jacobian and B
        % material[Material]: Requires methods E, nu, and D
        % order [Int]: Gauss integration order
            % Constitutive Relationship
            % Both material properties skip 3rd col because it is a plain
            % stress problem
            elastic = Physics.ElasticShell(material)
            piezo = material.D(:,[1 2 4 5 6]);  
            electric = diag(material.e);
            C = [   elastic piezo';
                    piezo   electric];
            % Function to be integrated
            function K_in_point = K_in_point(ksi,eta,zeta)
                % Piezo Part, generates the Electric Field (only z
                % component from 2 or 1 voltage element dof.
                jac = element.jacobian(ksi,eta,zeta);
                cosines = Element.direction_cosines(jac);
                inv_jac = jac \ eye(3);
                dof_per_ele = 1;
                if (dof_per_ele == 2)
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
                B = blkdiag(Element.T(jac)*B_mech,cosines*dN_xyz);
                K_in_point = B'*C*B*det(jac);
            end
            fun_in = @(xi,eta,mu) (K_in_point(xi,eta,mu));
            K = Integral.Volume3D(fun_in,order,[-1 1]);
        end
        function E_z = Electric_Field_z(element,ksi,eta,zeta)
        end
        function L = apply_volume_load(element,order,q)
        % L = apply_load(element,order,q)
        % Generates a load dof vector by integrating a constant load along
        % the element.
        % L [n_ele_dofs x 1][Float]: Load vector
        % element [Element]
        % order [Int]: Gauss integration order
        % q [dof x 1][Float]: Constant applied load
            function L_out = apply_load_in_point(ksi,eta,zeta)
                NN = Element.shape_to_diag(length(q),element.N(ksi,eta)); 
                L_out = NN'*det(element.jacobian(ksi,eta,zeta))*q;
            end
            fun_in = @(xi,eta,mu) (apply_load_in_point(xi,eta,mu));
            L = Integral.Volume3D(fun_in,order,[-1 1]);
        end
        function K = K_Shell(element,material,order)
        % K [n_dofxn_dof][Float] Stiffness as calculated in Cook 361 12.4-14
        % element [Element]: Requires methods jacobian and B
        % material[Material]: Requires methods E and nu
        % order [Int]: Gauss integration order
            C = Physics.ElasticShell(material);
            function K_in_point = K_in_point(ksi,eta,zeta)
                jac = element.jacobian(ksi,eta,zeta);
                B = Element.T(jac)* ...
                    Physics.B_Shell(element,ksi,eta,zeta); % Cook [7.3-10]
                K_in_point = B'*C*B*det(jac);
            end
            fun_in = @(xi,eta,mu) (K_in_point(xi,eta,mu));
            K = Integral.Volume3D(fun_in,order,[-1 1]);
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
            jac = element.jacobian(ksi,eta,zeta);
            invJac = jac \ eye(3);
            dN = invJac(:,1:2)*element.dN(ksi,eta);

            B  = zeros(6,element.n_nodes*dofs_per_node);
            % B matrix has the same structure for each node, written as
            % [aux1 aux2].
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
    