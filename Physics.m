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
        function K = K_Shell(element,material,order)
        % K [20x20][Float] Stiffness as calculated in Cook 361 12.4-14
        % element [Element]: Requires methods jacobian and B
        % material[Material]: Requires methods E and nu
        % order [Int]: Gauss integration order
            C = Physics.ElasticShell(material);
            function out = K_in_point(xi,eta,mu)
                % Compute Jacobian
                jacobian = element.jacobian(xi,eta,mu);
                % Linematic matrix for stiffness
                B = Physics.B(element,xi,eta,mu);            
                out = B'*C*B*det(jacobian);
                % POSIBLE ERROR: va con t multiplicando?
                % B'*D*B*wtx*wty*det(jacobian);
            end
            fun_in = @(xi,eta,mu) (K_in_point(xi,eta,mu));
            K = Integral.Volume3D(fun_in,20,order,[-1 1]);
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
            C(3,:) = 0;
            C(:,3) = 0;
            C(5,5) = (5/6)*C(5,5);
            C(6,6) = (5/6)*C(6,6);
        end
        function B_out = B(element,xi,eta,mu)
            % B_out [6x20] Computes B matix - Cook 361 12.5-10
            % B = H*T_j*N_devs
            B_out = Physics.H*Physics.invT(element,xi,eta,mu)*...
                                Physics.NdevsShell(element,xi,eta,mu);
        end     
        function H_out = H
            % H_out = H
            % H_out [6x9] Cook pg 181 6.7-5
            % goes from diff_U_xyz to Strain vector
            H_out = zeros(6,9);
            H_out([1 10 18 22 26 35 42 47 51]) = 1;
        end
        function T_out = invT(element,xi,eta,mu)
            % T_out [9x9][Float] Cook 360 12.5-8
            inv_jac = inv(element.jacobian(xi,eta,mu));
            T_out = blkdiag(inv_jac,inv_jac,inv_jac);
        end
        function Ndevs_out = NdevsShell(element,xi,eta,mu)
            % Ndevs_out = Ndevs(element,xi,eta,mu)
            % Ndevs_out [9x20] Cook pg 360 12.5-9 from Dofs U_i to diffU
            % It belongs in Physics because it contains the way the 
            % mechanical degrees of freedom combine to yield the 
            % displacement field
            dim = 3;
            ndofs = 5;
            n_nodes = element.n_nodes;
            t_at_node = element.thickness_at_node; 
            dN = [Element.dNdxi_Q4(xi,eta); Element.dNeta_Q4(xi,eta)];
            N = Element.N_Q4(xi,eta);
            mu_mat = element.mu_matrix;
            Ndevs_out = zeros(dim^2,n_nodes*ndofs);
            % The matrix in Cook is specified for one node with subscript i
            % The loop with 'node' computes each node matrix in aux, and then
            % stores it in Ndevs_out
            % The loop with 'j' computes aux every 3 rows, which correspond
            % to u, then v, and then w.
            for node = 1:n_nodes
                aux = zeros(dim^2,ndofs);
                for j = 1:dim
                    aux2 = zeros(dim,ndofs);
                    aux2(1:2,j) = dN(:,node);
                    aux2(1:2,4:5) = 0.5*t_at_node(node)*mu* ...
                                            dN(:,node)*mu_mat(j,:,node);
                    aux2(3,4:5) = 0.5*t_at_node(node)* ...
                                            N(node)*mu_mat(j,:,node);
                    aux(index_range(dim,j),:) = aux2;
                end 
                Ndevs_out(:,index_range(ndofs,node)) = aux;
            end
        end
    end
end
    