classdef Mesh
    properties
        coords
        connect
        normals
        thickness
        ele_type
    end
    properties (Dependent)
        n_nodes
        n_ele
        nodes_per_ele
    end
    methods (Static)
        function nodal_systems = nodal_systems(eleType,coords,connect)
            n_coords = size(coords,1);              % Number of coords
            ksi = [ -1  1  1 -1  0  1  0 -1  0 ];
            eta = [ -1 -1  1  1 -1  0  1  0  0 ];
            tol = 1e-9;
            nodal_systems = zeros(3,3,n_coords);
            for inod = 1:n_coords
                % ele:          Find the connect where the node is used
                % localNode:    Find which node it represents to the element
                [ele, localNode] = find(connect == inod);
                nuso = length(ele);     % Number of connect the node belongs to.
                % Loop through those connect to find v
                v = zeros(nuso,3);      % Vector directions (v3) of the node in each element
                for iele = 1:nuso
                    % Get the jacobian of the element
                    dN = Element.shapefunsder([ksi(localNode(iele)),eta(localNode(iele))],eleType);
                    elecoords = connect(ele(iele),:);
                    nodalCoords = coords(elecoords,:);
                    jac = dN*nodalCoords;
                    % Use the jacobian to get v3
                    v3 = cross(jac(1,:),jac(2,:));
                    v(iele,:) = v3/norm(v3);
                end
                v3 = mean(v,1);
                % from v3 get v1 and v2
                v3 = v3'/norm(v3);
                v1 = cross([0 1 0]',v3);
                largo = norm(v1);
                if largo > tol
                    v1 = v1/largo;
                else
                    v1 = [1 0 0]';
                end
                % v2 is orthogonal both to v1 and v3
                v2 = cross(v3,v1);
                v2 = v2/norm(v2);
                % save v1 v2 and v3 in the global coords' directions
                nodal_systems(:,:,inod) = [v1 v2 v3];
            end
        end
    end
    methods
        function obj = Mesh(ele_type,coords,connect,thickness_in)
            % obj = Mesh(coords,connect,node_normals)
            % Sanity Checks
            require(isnumeric(coords) && isnumeric(connect), ...
                'ArgumentError: all arguments should be numeric')
            require(all(all(~mod(connect,1))), ...
                'ArgumentError: connect should be all integers')
            require(max(max(connect))==size(coords,1), ...
                'ArgumentError: biggest node in connect should be size coords')
            % Setting variables
            obj.ele_type = ele_type;
            obj.coords = coords;
            obj.connect = connect;
            obj.normals = Mesh.nodal_systems(ele_type,coords,connect);
            obj.thickness = thickness_in;
        end
        function L = assembly_vector(mesh,dofs_per_node,dofs_per_ele,fun_in)
            % K = assembly_vector(mesh,dofs_per_node,dofs_per_ele,fun_in)
            % K [n_dofs x n_dofs]
            % dofs_per_node [Int]
            % dofs_per_ele [Int]
            % fun_in [Fhandle] f(element) -> [dofs_per_ele x 1]
            % fun_in follows convention [node_dofs ele_dofs]
            require(all(~mod([dofs_per_node,dofs_per_ele],1)), ...
                'ArgumentError: dofs should be integers');
            L = zeros(mesh.n_dofs(dofs_per_node,dofs_per_ele),1);
            % Loop through elements
            for e = 1:mesh.n_ele
                ele = mesh.ele(e);
                dofs = mesh.all_eles_dofs(dofs_per_node,dofs_per_ele,e);
                L(dofs) = L(dofs) + fun_in(ele);
            end
        end
        function K = assembly_matrix(mesh,dofs_per_node,dofs_per_ele,fun_in)
            % K = assembly_matrix(mesh,dofs_per_node,dofs_per_ele,fun_in)
            % K [n_dofs x n_dofs]
            % dofs_per_node [Int]
            % dofs_per_ele [Int]
            % fun_in [Fhandle] f(element) -> [dofs_per_ele x dofs_per_ele]
            % fun_in follows convention [node_dofs ele_dofs]
            require(all(~mod([dofs_per_node,dofs_per_ele],1)), ...
                'ArgumentError: dofs should be integers');
            K = zeros(mesh.n_dofs(dofs_per_node,dofs_per_ele));
            % Loop through elements
            for e = 1:mesh.n_ele
                ele = mesh.ele(e);
                dofs = mesh.all_eles_dofs(dofs_per_node,dofs_per_ele,e);
                K(dofs,dofs) = K(dofs,dofs) + fun_in(ele);
            end
        end
        function nodes = find_nodes(mesh,condition)
            % nodes = find_nodes(fun_in)
            % nodes [n_nodes x 1][Bool]
            % condition [FHandle][(x,y,z) -> Bool]
            % Returns a logical matrix where the true-valued nodes satisfy
            % the condition
            nodes = false(mesh.n_nodes,1);
            for n = 1:mesh.n_nodes
                x = mesh.coords(n,1);
                y = mesh.coords(n,2);
                z = mesh.coords(n,3);
                nodes(n) = condition(x,y,z);
            end
        end
        %% DOF helpers
        function dofs = all_eles_dofs(mesh,dofs_per_node,dofs_per_ele,ele_id)
            % eles_dofs(mesh,dofs_per_node,dofs_per_ele,ele_id)
            % Computes the complete Dof list for a certain element
            node_dofs = mesh.node_dofs(dofs_per_node,ele_id);
            ele_dofs = mesh.eles_dofs(dofs_per_node,dofs_per_ele,ele_id);
            dofs = [node_dofs; ele_dofs];
        end
        function dofs = eles_dofs(mesh,dofs_per_node,dofs_per_ele,ele_id)
            % dofs(mesh,dofs_per_node,dofs_per_ele,ele_id)
            % Finds the dofs that correspond to the element but not to the 
            % nodes of the element
            dofs = mesh.last_node_dof(dofs_per_node) + ...
                index_range(dofs_per_ele,ele_id)';
        end
        function dofs = node_dofs(mesh,dofs_per_node,ele_id)
            % dofs = node_dofs(mesh,dofs_per_node,ele_id)
            % Finds the dofs that correspond only to the nodes of the
            % element
            nodes = mesh.ele_nodes(ele_id);
            % Loop through nodes and add their dofs
            dofs = zeros(mesh.nodes_per_ele*dofs_per_node,1);
            for i = 1:mesh.nodes_per_ele
                dofs(index_range(dofs_per_node,i)) = ... 
                    index_range(dofs_per_node,nodes(i));
            end
        end
        function out = n_dofs(mesh,dofs_per_node,dofs_per_ele)
            % out = n_dofs(mesh,dofs_per_node,dofs_per_ele)
            % Total number of DOFs
            out = mesh.last_node_dof(dofs_per_node) + ...
                    mesh.n_ele * dofs_per_ele;
        end
        function out = last_node_dof(mesh,dofs_per_node)
            % out = last_node_dof(mesh,dofs_per_node)
            out = mesh.n_nodes*dofs_per_node;
        end
        %% Element Helpers
        function ele = ele(mesh,ele_id)
            % ele = ele(mesh,ele_id)
            % ele [Element]
            % Create element with ID ele_id from the mesh
            nodes = mesh.ele_nodes(ele_id);
%             v = mesh.normals(:,:,elecoords);     % Direction vectors from iele's coords
            ele = Element(mesh.ele_type,mesh.coords(nodes,:),mesh.normals(:,:,nodes), ...
                            mesh.thickness(nodes));
        end
        function node_ids = ele_nodes(mesh,ele_id)
            % nodes_ids = ele_nodes(mesh,ele_id)
            % nodes_ids [1xnodes_per_element]
            % Returns all the IDs corresponding to ele_id
            node_ids = mesh.connect(ele_id,:);
        end
        %% Dependent Properties
        function out = get.n_nodes(mesh)
        % out = n_nodes(mesh)
        % Number of Nodes
            out = size(mesh.coords,1);
        end
        function out = get.n_ele(mesh)
        % out = n_nodes(mesh)
        % Number of Elements
            out = size(mesh.connect,1);
        end
        function out = get.nodes_per_ele(mesh)
        % out = nodes_per_ele(mesh)
            out = size(mesh.connect,2);
        end
    end
end