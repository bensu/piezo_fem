classdef Element
    properties
        coords
        normals
    end
    properties (Dependent)
        n_nodes
        thickness_at_node
        mu_matrix
    end
    methods
        function obj = Element(coords,normals)
            require(size(coords,1)==4, ...
                'ArgumentError: only 4 nodes');
            require(size(coords)==size(normals), ...
                'ArgumentError: coords and normals should have same size');
            obj.coords = coords;
            obj.normals = normals;
        end
        function jac_out = jacobian(element,xi,eta,mu)
            % jac_out = jacobian(element,xi,eta,mu)
            % jac_out [3x3][Float]: Jacobian Matrix
            % xi, eta, mu [Float] between [-1,1]
            % Computes the Jacobian for the ShellQ4 element
            % as defined in Cook pg 360 12.5-4
            dNdxi = Element.dNdxi_Q4(xi,eta);
            dNdeta = Element.dNeta_Q4(xi,eta);
            N = Element.N_Q4(xi,eta);
            % This ASSUMES element.normals has a certain shape
            jac_out = [dNdxi;dNdeta;zeros(1,4)]*element.coords ...
                + [mu*dNdxi;mu*dNdeta;N]*element.normals/2;
        end
        function thickness_out = get.thickness_at_node(element)
            % thickness_out = get.thickness_at_node(element)
            % thickness_out [4x1][Float] original thickness of the shell at
            % each node. Cook pg 358 12.5-1 as t_i
            % Computed from the normals vector
            V3 = element.normals';
            node_num = 4;
            thickness_out = zeros(node_num,1);
            for node = 1:node_num
                thickness_out(node) = norm(V3(:,node));
            end
        end
        function mu_out = get.mu_matrix(element)
            % mu_out = get.mu_matrix(element)
            % mu_out [3x2x4] Cook pg 359 12.5-3
            % Computes the matrix for each element and stores it in a 3D
            % array
            node_num = 4;
            V1 = zeros(3,node_num);
            node = 1;
            for i = [-1,1]
                for j = [-1,1]
                    V1(:,node) = Element.dNeta_Q4(i,j)*element.coords;
                    node = node + 1;
                end
            end
            V3 = element.normals';
            mu_out = zeros(3,2,node_num);
            for node = 1:node_num
                mu_out(:,2,node) = V1(:,node)/norm(V1(:,node)); 
                aux = cross(V3(:,node),mu_out(:,2,node));
                mu_out(:,1,node) = -aux/norm(aux);
            end
        end
        function out = get.n_nodes(element)
            % out = get.n_nodes(element)
            % Number of nodes in the element
            out = size(element.coords,1);
        end
    end
    methods (Static)
        function N_out = N_ShellQ4(element,xi,eta,mu)
            % Not really used!!!
            require(isnumeric(mu), ...
                'ArgumentError: xi, eta, and mu should be numeric')
            require(-1<=mu && mu<=1, ... 
                'ArgumentError: mu should is not -1<=mu<=1')
            % Need to check the way it works with Jacobian
            N_out = Element.N_Q4(xi,eta)*mu*element.normals;
        end
        function N_out = N_Q4(xi,eta)
            %  Notes:
            %     1st node at (-1,-1), 3rd node at (-1,1) 
            %     4th node at (1,1), 2nd node at (1,-1)
            require(isnumeric([xi eta]), ...
                'ArgumentError: Both xi and eta should be numeric')
            require(-1<=xi && xi<=1, ...
                'ArgumetnError: xi should be -1<=xi<=1')
            require(-1<=eta && eta<=1, ...
                'ArgumetnError: eta should be -1<=eta<=1')
            N_out = zeros(1,4);
            N_out(1) = 0.25*(1-xi)*(1-eta);
            N_out(3) = 0.25*(1+xi)*(1-eta);
            N_out(2) = 0.25*(1-xi)*(1+eta);
            N_out(4) = 0.25*(1+xi)*(1+eta);
        end
        function dNdxiQ4_out = dNdxi_Q4(xi,eta)
            % derivatives
            require(isnumeric([xi eta]), ...
                'ArgumentError: Both xi and eta should be numeric')
            require(-1<=xi && xi<=1, ...
                'ArgumetnError: xi should be -1<=xi<=1')
            require(-1<=eta && eta<=1, ...
                'ArgumetnError: eta should be -1<=eta<=1')
            dNdxiQ4_out = zeros(1,4);
            dNdxiQ4_out(1) = -0.25*(1-eta);
            dNdxiQ4_out(2) = 0.25*(1-eta);
            dNdxiQ4_out(3) = -0.25*(1+eta);
            dNdxiQ4_out(4) = 0.25*(1+eta);
        end
        function dNdetaQ4_out = dNeta_Q4(xi,eta)
            % derivatives
            require(isnumeric([xi eta]), ...
                'ArgumentError: Both xi and eta should be numeric')
            require(-1<=xi && xi<=1, ...
                'ArgumetnError: xi should be -1<=xi<=1')
            require(-1<=eta && eta<=1, ...
                'ArgumetnError: eta should be -1<=eta<=1')
            dNdetaQ4_out = zeros(1,4);
            dNdetaQ4_out(1) = -0.25*(1-xi);
            dNdetaQ4_out(2) = -0.25*(1+xi);
            dNdetaQ4_out(3) = 0.25*(1-xi);
            dNdetaQ4_out(4) = 0.25*(1+xi);
        end
    end
end
            