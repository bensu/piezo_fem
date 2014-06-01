classdef Factory
    % Factory
    % Use to generate common cases for tests, for example meshes and
    % elements that can be later used.It can provide parameters and memory 
    % environments when needed.
    % In principle, all methods should be static, since the class works as a
    % library and not as a memory structure.
    methods (Static)
        function mesh = BrickMesh(type,n_elements,sides)
            % mesh = BrickMesh(type,n_elements,sides)
            % mesh [Mesh]: New Generated mesh
            % sides = [a,b,c] [1x3][Float]: sides of the shell 
            % n_elements = [m,n,o] [1x2]3Int]: num of elements in each edge
            [coords, connect] = Factory.Box(n_elements,sides);
            thickness = sides(3)*ones(1,size(coords,1));
            mesh = Mesh(type,coords,connect,thickness);       
        end
        function mesh = ShellMesh(type,n_elements,sides)
            % mesh = ShellMesh(n_elements,sides)
            % mesh [Mesh]: New Generated mesh
            % sides = [a,b,c] [1x3][Float]: sides of the shell 
            % n_elements = [m,n] [1x2][Int]: num of elements in each edge
            [coords, connect] = Factory.Plate(type,n_elements,sides);
            thickness = sides(3)*ones(1,size(coords,1));
            mesh = Mesh(type,coords,connect,thickness);
        end
        function [coords, connect, node_normals] = Shell(n_elements,sides)
            % mesh = ShellMesh(n_elements,sides,E,nu,rho)
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
        function [coords, connect] = Plate(type,n_elements,sides)
            % [coords, connections] = Plate(n_elements,sides)
            % coords [nodes x 3][Float]: coords of all the nodes
            % connections [elements x node_per_ele][Int]: 
            %       node_ids associated to each element
            % n_elements [1x3][Int]: number of elements along each
            % dimension
            % sides [1x3][Float]: length of plate box in each dimension
            %   side(3) = z : where is the plate in space.
            % Sanity Checks
            require(all(~rem(n_elements,1)), ...
                'ArgumentError: n_elements should be integers')
            require(isnumeric(sides), ...
                'ArgumentError: sides should be numeric')
            require(length(sides)==3 && length(n_elements)==2, ...
                'ArgumentError: sides, n_element should be [3x1] & [2x1]');
            % Prepare Variables
            m = n_elements(1);  % Number of elements in x direction
            n = n_elements(2);  % Number of elements in y direction
            nx = 2*m + 1; % Number of nodes in x direction
            ny = 2*n + 1; % Number of nodes in y direction
            dx = sides(1)/(nx-1);    % Distance between points in x
            dy = sides(2)/(ny-1);    % Distance between points in y
            z = 0;            % Height where the plate is placed
            
            % Generate coordinates
            coords = zeros(nx*ny,3);
            n_count = 1;
            for j = 1:ny
                for i = 1:nx
                    coords(n_count,:) = [(i-1)*dx,(j-1)*dy,z];
                    n_count = n_count + 1;
                end
            end
            
            % Generate connectivity
            a = 1:3;    % Helper variable
            connect = zeros(m*n,9);
            e_count = 1;    % Element count
            for j = 1:n
                for i = 1:m
                    connect(e_count,:) = [2*(i-1)+2*nx*(j-1)        + a ...
                                          nx+2*(i-1)+2*nx*(j-1)     + a ...
                                          2*nx+2*(i-1)+2*nx*(j-1)   + a];
                    e_count = e_count + 1;
                end
            end
            connect = connect(:,[1 3 9 7 2 6 8 4 5]);
            switch type
                case {EleType.Q4 EleType.AHMAD4}
                    [coords, connect] = Factory.remove_nodes(coords,connect(:,1:4));
                case {EleType.Q8 EleType.AHMAD8}
                    [coords, connect] = Factory.remove_nodes(coords,connect(:,1:8));
                case {EleType.Q9 EleType.AHMAD9}
                otherwise
                    error('Wrong type of element: should be Q4, Q8, or Q9');
            end
        end
        function [coords, connect] = remove_nodes(coords,connect)
        % [coords, connect] = remove_nodes(coords,connect)
        % Coords leaves only the nodes present in connect, and then renames
        % those of connnect with their new order in coords
            present_nodes = unique(reshape(connect,[],1));
            coords = coords(present_nodes,:);
            new_connect = zeros(size(connect));
            for i = 1:length(present_nodes)
                new_connect(connect == present_nodes(i)) = i;
            end
            connect = new_connect;
        end
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
            require(length(sides)==3 && length(n_elements)==3, ...
                'ArgumentError: sides and n_element should be [3x1]');
            
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
            node_count = 1; % Node count
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
            e = 1;  % Element count
            for k = 1:p
                for j = 1:n
                    for i = 1:m
                        connections(e,1) = 1+(i-1)+(j-1)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections(e,2) = 1+i+(j-1)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections(e,3) = i+(j)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections(e,4) = 1+i+(j)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections(e,5) = 1+(i-1)+(j-1)*(m+1)+k*(m+1)*(n+1);
                        connections(e,6) = 1+i+(j-1)*(m+1)+k*(m+1)*(n+1);
                        connections(e,7) = i+(j)*(m+1)+k*(m+1)*(n+1);
                        connections(e,8) = 1+i+(j)*(m+1)+k*(m+1)*(n+1);
                        e = e + 1;
                    end
                end
            end
           connections = connections(:,[1 2 4 3 5 6 8 7]);
        end
    end
end