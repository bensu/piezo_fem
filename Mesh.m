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
            obj.normals = obj.nodal_systems();
            obj.thickness = thickness_in;
        end
        function nodal_systems = nodal_systems(mesh)
            ksi = [ -1  1  1 -1  0  1  0 -1  0 ];
            eta = [ -1 -1  1  1 -1  0  1  0  0 ];
            tol = 1e-9;
            nodal_systems = zeros(3,3,mesh.n_nodes);
            for n = 1:mesh.n_nodes
                % ele:          Find the connect where the node is used
                % localNode:    Find which node it represents to the element
                [ele, localNode] = find(mesh.connect == n);
                nuso = length(ele);     % Number of connect the node belongs to.
                % Loop through those connect to find v
                v = zeros(nuso,3);      % Vector directions (v3) of the node in each element
                element = Element(1,mesh.ele_type,[],[],[]);
                for e = 1:nuso
                    % Create the element
                    ele_nodes = mesh.connect(ele(e),:);
                    ele_coords = mesh.coords(ele_nodes,:);
                    
                    % Get some simplified Jacobian of the element
                    dN = element.dN(ksi(localNode(e)),eta(localNode(e)));
                    jac = dN*ele_coords;
                    
                    % Use the jacobian to get v3
                    v3 = cross(jac(1,:),jac(2,:));
                    v(e,:) = v3/norm(v3);
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
                nodal_systems(:,:,n) = [v1 v2 v3];
            end
        end
        function L = integral_along_surface(mesh,dofs_per_node, ...
                                        dofs_per_ele,s_nodes,fun_in)
        % L = integral_along_surface(mesh,s_nodes,fun_in)
        % s_nodes [n x 1]: nodes in the surface
        % fun_in [FHandle](element,surface_c,surface_value)
            ele_ids = mesh.unique_eles(s_nodes);
            ele_s = mesh.ele_surfaces(s_nodes,ele_ids);
            [~, s_coord, s_val] = Element.surfaces(mesh.ele_type);
            function L_out = aux_fun(element)
                surfaces = ele_s{element.id == ele_ids};    % Element's surface #
                L_out = [];
                for s = surfaces
                    if isempty(L_out)
                        L_out = fun_in(element,s_coord(s),s_val(s));
                    else
                        L_out = L_out + fun_in(element,s_coord(s),s_val(s));
                    end
                end
            end
            ele_fun = @(element) aux_fun(element);
            L = mesh.assembly_vector_along(dofs_per_node,dofs_per_ele, ...
                                                        ele_ids,ele_fun);
        end
        function L = assembly_vector_along(mesh,dofs_per_node,dofs_per_ele, ...
                                                            ele_ids,fun_in)
        % L = assembly_vector_along(mesh,dofs_per_node,dofs_per_ele, ...
        %                                    ele_ids,fun_in)
        % K [n_dofs x n_dofs]
        % dofs_per_node [Int]
        % dofs_per_ele [Int]
        % fun_in [Fhandle] f(element) -> [dofs_per_ele x 1]
        % fun_in follows convention [node_dofs ele_dofs]
            require(all(~mod([dofs_per_node,dofs_per_ele],1)), ...
                'ArgumentError: dofs should be integers');
            mesh.n_dofs(dofs_per_node,dofs_per_ele)
            L = sparse(mesh.n_dofs(dofs_per_node,dofs_per_ele),1);
            % Loop through elements
            for e = ele_ids
                ele = mesh.ele(e);
                dofs = mesh.all_eles_dofs(dofs_per_node,dofs_per_ele,e);
                L(dofs) = L(dofs) + fun_in(ele);
            end
        end
        function L = assembly_vector(mesh,dofs_per_node,dofs_per_ele,fun_in)
        % L = assembly_vector(mesh,dofs_per_node,dofs_per_ele,fun_in)
        % L [n_dofs x 1]
        % dofs_per_node [Int]
        % dofs_per_ele [Int]
        % fun_in [Fhandle] f(element) -> [dofs_per_ele x 1]
        % fun_in follows convention [node_dofs ele_dofs]
            L = assembly_vector_along(mesh,dofs_per_node,dofs_per_ele, ...
                                                    1:mesh.n_ele,fun_in);
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
            n_dofs = mesh.n_dofs(dofs_per_node,dofs_per_ele);
            K = sparse(n_dofs,n_dofs);
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
                r = mesh.coords(n,:);   % r = [x,y,z]
                nodes(n) = condition(r(1),r(2),r(3));
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
            ele = Element(ele_id,mesh.ele_type,mesh.coords(nodes,:), ...
                        mesh.normals(:,:,nodes),mesh.thickness(nodes));
        end
        function node_ids = ele_nodes(mesh,ele_id)
        % nodes_ids = ele_nodes(mesh,ele_id)
        % nodes_ids [1xnodes_per_element]
        % Returns all the IDs corresponding to ele_id
            node_ids = mesh.connect(ele_id,:);
        end
        function ele_ids = node_eles(mesh,node_id)
        % ele_ids = node_eles(mesh,node_id)
        % node_id [n_nodes x 1][Int]: Ids of the nodes.
        % ele_ids {n_nodes x 1}[Int]: Ids of the elements that contain
        % those nodes.
            % A cell array is used because each node might return a
            % different number of elements
            n_nodes = length(node_id);
            ele_ids = cell(n_nodes,1);
            for n = 1:n_nodes
                ele_ids{n} = find(any((mesh.connect == node_id(n))')')';
            end
        end
        function ele_ids = unique_eles(mesh,node_ids)
            % ele_ids = unique_elements(mesh,node_ids)
            % node_id [n_nodes x 1][Int]: Ids of the nodes.
            % ele_ids [? x 1][Int]: Ids of the elements that contain those nodes.
            node_eles = mesh.node_eles(node_ids);
            ele_ids = [];
            for n = 1:length(node_eles)
                ele_ids = [ele_ids node_eles{n}];
            end
            ele_ids = unique(ele_ids);
        end
        function ele_surfaces = ele_surfaces(mesh,nodes,ele_ids)
            % ele_surfaces = ele_surfaces(mesh,nodes,ele_ids)
            % ele_surfaces {n_ele_ids x 1}[Int] the surfaces (between 1
            % and 6) of each element that contain the nodes.
            % nodes: [n x 1][Int] node_ids of the nodes in the global surface
            n_ele = length(ele_ids);    % Number of elements
            ele_surfaces = cell(n_ele,1);
            % surfaces contains data from the element type
            surfaces = Element.surfaces(mesh.ele_type);
            % Loop through all elements
            for e = 1:n_ele
                % Get all of its nodes
                ele_nodes = mesh.ele_nodes(ele_ids(e));
                % Loop through the nodes
                ele_surface_nodes = false(1,length(ele_nodes));
                for n = 1:length(ele_nodes)
                    % Check if the node is in the global surface
                    ele_surface_nodes(n) = any(nodes == ele_nodes(n));
                end
                % Those who were, are stored in ele_surface_nodes
                ele_surface_nodes = sort(find(ele_surface_nodes));
                ele_surface = [];
                % Check which surfaces contain ele_surface_nodes
                for s = 1:size(surfaces,2)
                    if all(ele_surface_nodes == surfaces(s,:))
                        ele_surface = [ele_surface s];
                    end
                end
                ele_surfaces{e} = ele_surface;
            end
        end
        %% Plotting
        function plot(mesh)
            % Nodes per element (to be plotted)
            nodes_per_ele = 4;  % HARDCODED FOR RECTANGULAR SHELL ELEMENTS                              
            % Initialization of the required matrices
            X = zeros(nodes_per_ele,mesh.n_ele); 
            Y = zeros(nodes_per_ele,mesh.n_ele); 
            Z = zeros(nodes_per_ele,mesh.n_ele);
            for e = 1:mesh.n_ele
                e_nodes = mesh.ele_nodes(e);    % nodes for ele with id = e
                for n = 1:nodes_per_ele
                    X(n,e) = mesh.coords(e_nodes(n),1);    % extract x value of the node
                    Y(n,e) = mesh.coords(e_nodes(n),2);    % extract y value of the node
                    Z(n,e) = mesh.coords(e_nodes(n),3);
                end
            end
            
            % Plotting the FEM mesh, diaplay Node numbers and Element numbers
            f1 = figure ;
            set(f1,'name','Mesh','numbertitle','off') ;
            patch(X,Y,Z,'w');
            hold on
            e_nodes = mesh.connect(:,1:end)';
            for n = 1:mesh.n_ele
                text(X(:,n),Y(:,n),Z(:,n),int2str(e_nodes(1:4,n)), ...
                                                'fontsize',8,'color','k');
                text(sum(X(:,n))/4,sum(Y(:,n))/4,sum(Z(:,n))/4,int2str(n), ...
                                                'fontsize',10,'color','r') ;
            end
            hold off
        end
        function plot_displacements(mesh,scale,displacements)
            % plot_displacements(mesh,scale,displacements)
            % displacements [n_nodes x 3]
            % scale [Float]: the bigger scale, the bigger displacements look
            aux_mesh = mesh;
            aux_mesh.coords = mesh.coords + scale*displacements;
            aux_mesh.plot();
        end
        function plot_vector_field(mesh,vector_field,options)
            % plot_vector_field(mesh,vector_field,option)  
            % Wrapper for quiver3
            % vector_field [n_nodes x 3][Float]
            % options: MATLAB plotting options, i.e. 'k-'
            quiver3(mesh.coords(:,1),mesh.coords(:,2),mesh.coords(:,3), ...
            	vector_field(:,1),vector_field(:,2),vector_field(:,3), ...
                options);
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