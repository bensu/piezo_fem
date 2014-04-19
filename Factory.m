classdef Factory
    methods (Static)
        function [coords, connections] = Box(n_elements,sides)
            % [coords, connections] = BasicBox(n_elements,sides)
            % coords [nodes x 3][Float]: coords of all the nodes
            % connections [elements x 8][Int]: node_ids associated to each element
            % n_elements [1x3][Int]: number of elements along each
            % dimension
            % sides [1x3][Float]: length of the box in each dimension
            
            % Sanity Checks
            require(all(~rem(n_elements,1)), ...
                'ArgumentError: n_elements should be integers')
            require(isnumeric(sides), ...
                'ArgumentError: sides should be numeric')
            
            % Prepare Variables
            m = n_elements(1);  % Number of elements in x direction
            n = n_elements(2);  % Number of elements in y direction
            p = n_elements(3);  % Number of elements in z direction
            % Distance between nodes in each dimension
            dx = sides(1)/m;
            dy = sides(2)/n;
            dz = sides(3)/p;
            coords = zeros((n+1)*(m+1)*(p+1),3);
            connections = zeros(n*m*p,8);
            
            % Create all coords by looping
            % We loop along the x,y,z dimensions, yet the data structure
            % coords grows linearly its rows.
            % node_count keeps track of the row to be updated.
            node_count = 1;
            for k = 1:p+1
                for j = 1:n+1
                    for i = 1:m+1
                        coords(node_count,:) = [(i-1)*dx,(j-1)*dy,(k-1)*dz];
                        node_count = node_count + 1;
                    end
                end
            end
            
            % Create all connections by looping,
            % works similarly to coords
            ele_count = 1;
            for k = 1:p
                for j = 1:n
                    for i = 1:m
                        connections(ele_count,1) = 1+(i-1)+(j-1)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections(ele_count,2) = 1+i+(j-1)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections(ele_count,3) = i+(j)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections(ele_count,4) = 1+i+(j)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections(ele_count,5) = 1+(i-1)+(j-1)*(m+1)+k*(m+1)*(n+1);
                        connections(ele_count,6) = 1+i+(j-1)*(m+1)+k*(m+1)*(n+1);
                        connections(ele_count,7) = i+(j)*(m+1)+k*(m+1)*(n+1);
                        connections(ele_count,8) = 1+i+(j)*(m+1)+k*(m+1)*(n+1);
                        ele_count = ele_count + 1;
                    end
                end
            end
        end
        function [coords, connect, node_normals] = Shell(n_elements,sides)
            % mesh = ShellMesh(n_elements,sides,E,nu,rho)
            % mesh [Mesh]: New generated mesh
            % sides = [a,b,c] [1x3][Float]: sides of the shell 
            % n_elements = [m,n] [1x2][Int]: num of elements in each edge
            
            % Create Main Data
            % Since BasicBox generates a 3D box and we need a Shell, we
            % specify only 1 element in the z direction.
            [coords_aux, connect_aux] = Factory.Box([n_elements,1],sides);

            % Generate Mesh from Box
            n_original_nodes = size(coords_aux,1);
            require(rem(n_original_nodes,2)==0,'Uneven number of nodes')
            
            % Find bottom an top nodes coords
            % (I'm ASSUMING they have a specific order, the last ones are
            % the bottom ones, depends on BasicBox)
            r_bottom = coords_aux(1:n_original_nodes/2,:); 
            r_top = coords_aux((n_original_nodes/2+1):end,:);
            coords = (r_bottom + r_top)/2;    % Middle coords
            
            % Because we reduced the node quantity we don't need all the 
            % nodes in connection. I ASSUME that the relevant ones are
            % located in 1:4 (Depends on BasicBox)
            connect = connect_aux(:,1:4);
            
            % Vector that gives the orientation of the shell for each node
            node_normals = (r_top - r_bottom); % [nodes x 3]
        end
    end
end