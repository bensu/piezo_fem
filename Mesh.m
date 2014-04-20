classdef Mesh
    properties
        coords
        connect
        normals
    end
    properties (Dependent)
        n_nodes
        n_ele
        nodes_per_ele
    end
    methods
        function obj = Mesh(coords,connect,node_normals)
            % obj = Mesh(coords,connect,node_normals)
            % Sanity Checks
            require(size(coords)==size(node_normals), ...
                'ArgumentError: coords and normals should be consistent')
            require(isnumeric(coords) && isnumeric(connect) &&  ... 
                isnumeric(node_normals), ...
                'ArgumentError: all arguments should be numeric')
            require(all(all(~mod(connect,1))), ...
                'ArgumentError: connect should be all integers')
            require(max(max(connect))==size(coords,1), ...
                'ArgumentError: biggest node in connect should be size coords')
            % Setting variables
            obj.coords = coords;
            obj.connect = connect;
            obj.normals = node_normals;
        end
        function K = assembly(mesh,dofs_per_node,dofs_per_ele,fun_in)
            % K = assembly(mesh,dofs_per_node,dofs_per_ele,fun_in)
            % K [n_dofs x n_dofs]
            % dofs_per_node [Int]
            % dofs_per_ele [Int]
            % fun_in [Fhandle] f(element) -> [dofs_per_ele x dofs_per_ele]
            % fun_in follows convention [node_dofs ele_dofs]
            require(all(~mod([dofs_per_node,dofs_per_ele],1)), ...
                'ArgumentError: dofs should be integers');
            K = zeros(mesh.n_dofs(dofs_per_node,dofs_per_ele));
            for e = 1:mesh.n_ele
                ele = mesh.ele(e);
                dofs = mesh.all_eles_dofs(dofs_per_node,dofs_per_ele,e);
                K(dofs,dofs) = K(dofs,dofs) + fun_in(ele);
            end
        end
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
        function ele = ele(mesh,ele_id)
            % ele = ele(mesh,ele_id)
            % ele [Element]
            % Create element with ID ele_id from the mesh
            nodes = mesh.ele_nodes(ele_id);
            ele = Element(mesh.coords(nodes,:),mesh.normals(nodes,:));
        end
        function node_ids = ele_nodes(mesh,ele_id)
            % nodes_ids = ele_nodes(mesh,ele_id)
            % nodes_ids [1xnodes_per_element]
            % Returns all the IDs corresponding to ele_id
            node_ids = mesh.connect(ele_id,:);
        end
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