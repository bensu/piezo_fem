classdef Element
    % Element Class
    % Works as a general data structure for an element. In principle, it
    % should not contain methods specific to a certain element type, but
    % raher the general ones.
    properties
        id          % [Int]: Unique id number sed by the mesh to represent the ele.
        type        % [EleType]: Contains all the specific functions.
        coords      % [n_nodes x 3][Float]: Node coordinates, in order
        normals     % [n_nodes x 3][Float]: Node coordinate system - Shell only
        thickness   % [n_nodes x 1][Float]: Thickness at each node - Shell only
        laminate    % [Laminate]
    end
    properties (Dependent)
        n_nodes     % [Int]: Number of nodes
        v3          % [3 x n_nodes]: Vector perpendicular to the shell at nodes - Shell only
    end
    methods
        function obj = Element(id,type,laminate,coords,normals,t_in)
            % function obj = Element(type,coords,normals,t_in)
            % Creates an element object
            %             require(size(coords,1)==4, ...
            %                 'ArgumentError: only 4 nodes');
            %             require(size(coords)==size(normals), ...
            %                 'ArgumentError: coords and normals should have same size');
            obj.id = id;
            obj.coords = coords;
            obj.thickness = t_in;
            obj.normals = normals;
            obj.type = type;
            obj.laminate = laminate;
        end
        function jac = jacobian(element,ksi,eta,zeta)
            % jac_out = jacobian(element,ksi,eta,zeta)
            % jac_out [3x3][Float]: Jacobian Matrix
            % ksi, eta, zeta [Float] between [-1,1], checking done in N
            % Computes the jacobian for Shell Elements
            % Cook [6.7-2] gives Isoparametric Jacobian
            % Cook [12.5-4] & [12.5-2] gives Shells Derivatives.
            switch element.type
                case {EleType.AHMAD4, EleType.AHMAD8, EleType.AHMAD9}
                    jac = element.ShellJac(ksi,eta,zeta);
                case {EleType.Q4,EleType.Q8,EleType.Q9}
                    jac = element.Jac2D(ksi,eta);
                case {EleType.H8}
                    jac = element.Jac3D(ksi,eta,zeta);
            end
        end
        function jac = ShellJac(element,ksi,eta,zeta)
            N  = element.N(ksi,eta,zeta);
            dN = element.dN(ksi,eta,zeta);
            t = element.thickness;
            tt = [t; t; t];
            v3t = (element.v3.*tt)';
            jac = [ dN*(element.coords + zeta*v3t/2);
                N*(v3t)/2 ];
        end
        function jac = Jac2D(element,ksi,eta)
            jac = element.dN(ksi,eta,0)*element.coords;
        end
        function jac = Jac3D(element,ksi,eta,zeta)
            jac = element.dN(ksi,eta,zeta)*element.coords;
        end
        function N = N(element,ksi,eta,zeta)
            % N = N(element,ksi,eta)
            % Element Shape Functions
            % Works by fetching them from EleType
            
            % Maybe a varargin should be implemented to avoid passing zeta
            % when element type is 2D.
            if Element.parameters_check(ksi,eta,zeta)
                N = element.type.N(ksi,eta,zeta);
            end
        end
        function dN = dN(element,ksi,eta,zeta)
            % N = dN(element,ksi,eta)
            % Element Shape Functions Derivatives
            % Works by fetching them from EleType
            if Element.parameters_check(ksi,eta,zeta)
                dN = element.type.dN(ksi,eta,zeta);
            end
        end
        function mu = mu_matrix(element)
            % mu = mu_matrix(cosines)
            % [mu] [3x2xn_nodes][Float] as defined in Cook 12.5-3
            % Specific Shell
            V = element.normals;
            V(:,3,:) = [];
            V = V(:,[2 1],:);
            V(:,1,:) = -V(:,1,:);
            mu = V;
        end
        function N = ShellN(element,ksi,eta,zeta)
            % N = N3(element,ksi,eta,zeta)
            % Shell proper shape functions
            % Specific Shell
            N = zeros(3,5*element.n_nodes);
            I = eye(3);
            mu = element.mu_matrix;
            t  = element.thickness;
            N2 = element.N(ksi,eta,0);
            for n = 1:element.n_nodes
                N(:,index_range(5,n)) = N2(n)*[I 0.5*t(n)*zeta*mu(:,:,n)];
            end
        end
        %% Dependent properties
       	function out = get.n_nodes(element)
            % out = get.n_nodes(element)
            % Number of nodes in the element
            out = size(element.coords,1);
        end
        function out = get.v3(element)
            out = squeeze(element.normals(:,3,:));
        end
    end
    methods (Static)
        function cosines = direction_cosines(jac)
            % local coordinate system [ksi eta zeta]
            dir1 = jac(1,:);
            dir3 = cross(dir1,jac(2,:));
            dir2 = cross(dir3,dir1);
            cosines = [ dir1/norm(dir1); dir2/norm(dir2); dir3/norm(dir3)];
        end
        function T = T(cosines)
            % jac [3x3][Float]: Jacobian matrix
            % Transformation of Strain, Cook pg 212: 
            % Cook [7.3-5]
            M1 = cosines;
            M2 = M1(:,[2 3 1]);
            M3 = M1([2 3 1],:);
            M4 = M2([2 3 1],:);
            T = [ M1.^2     M1.*M2;
                  2*M1.*M3  M1.*M4 + M3.*M2 ];
            % Since sigma_zz is ignored, we eliminate the appropriate row.
            T(3,:) = [];
        end
        function NN = shape_to_diag(dim,N)
            % NN = shape_to_diag(dim,N)
            % NN [Float][dim x n_nodes]: Repeated values of N
            % N [Float] [1 x n_nodes]: Evaluated shape functions
            % dim [Int]: dimension of the problem, 2 or 3.
            % Rearranges for some surface integral for loads
            n_nodes = length(N);
            NN = zeros(dim,n_nodes);
            I = eye(dim);
            for n = 1:n_nodes
                i = index_range(dim,n);
                NN(:,i) = N(n)*I;
            end
        end
        function bool = parameters_check(ksi,eta,zeta)
            % Throws error if any paramter is off
            % bool has no use.
            bool = false;
            require(isnumeric([ksi eta zeta]), ...
                'ArgumentError: Both ksi and eta should be numeric')
            require(-1<=ksi && ksi<=1, ...
                'ArgumetnError: ksi should be -1<=ksi<=1')
            require(-1<=eta && eta<=1, ...
                'ArgumetnError: eta should be -1<=eta<=1')
            require(-1<=zeta && zeta<=1, ...
                'ArgumetnError: zeta should be -1<=eta<=1')
            bool = true;
        end

    end
end
