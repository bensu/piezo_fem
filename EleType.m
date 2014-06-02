classdef EleType
    properties
        n_nodes             % nodes_per_element
        nodes_per_surface
    end
    enumeration
        Q4      (4,2)
        Q8      (8,3)
        Q9      (9,3)
        AHMAD4  (4,2)
        AHMAD8  (8,3)
        AHMAD9  (9,3)
        H8      (8,4)
    end
    methods
        function et = EleType(n_nodes,nodes_per_surface)
            % et = EleType(n_nodes,nodes_per_surface)
            % Initializer, should contain all the properties, for
            % enumeration
            et.n_nodes = n_nodes;
            et.nodes_per_surface = nodes_per_surface;
        end
        function [surfaces, s_coord, s_values] = surfaces(ele_type)
            % [surfaces, s_directions, s_values] = surfaces(ele_type)
            % Returns the nodes in each surface for the element type
            % Each element has 4 surfaces, each one has an id (1 through 4)
            % Each surface has either 1 or 3 nodes depending on the element.
            % The surface can be defined by fixing one of the local
            % coordinates, i.e. eta = 1 or ksi = -1.
            % surfaces [4 x 3][Int]: each row has the nodes of the surface
            % s_coord [4 x 1][Int]: each row has the local coordinate
            % that remains fixed to represent that surface.
            % s_values [4-6 x 1][Float]: the value of the coord that remains
            % fixed, usually -1 or 1.
            
            % BREAKS FOR H8
            surfaces = [1 2 5; 
                        2 3 6;
                        3 4 7;
                        1 4 8];
            s_coord     = [2 1 2 1]';
            s_values    = [-1 1 1 -1]';
            surfaces = surfaces(:,1:ele_type.nodes_per_surface);
        end
    end
    methods (Static)
        %% Shell Elements
        function types = Shell_Types()
           % types = Shell_Types()
           % types [1xn][EleType]: List with all the shell element types 
           types = [EleType.AHMAD4, EleType.AHMAD8, EleType.AHMAD9];
        end
        function bool = is_shell(ele_type)
            % bool = is_shell(ele_type)
            % Returns true if the ele_type argument is a shell.
            shell_types = EleType.Shell_Types();
            bool = false;
            for i = 1:length(shell_types)
                bool = bool || ele_type == shell_types(i);
            end
        end
        %% 3D elements
        %% H8 Methods
        function N = N_H8(ksi,eta,zeta)
            % N = N_H8(ksi,eta,zeta)
            % N [1 x 8][Float]: Shape Functions for H8 element
            % Follows anti-clockwise node-numbering convention:
            % (ksi,eta,zeta) = [(-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1)
            %                   (-1,-1,1),  (1,-1,1),  (1,1,1),  (-1,1,1)]
            N = zeros(1,8);
            count = 1;
            for k = [-1 1]
                for j = [-1 1]
                    for i = [-1 1]
                        N(count) = (1 + i*ksi)*(1 + j*eta)*(1 + k*zeta);
                        count = count + 1;
                    end
                end
            end
            N = N(:,[1 2 4 3 5 6 8 7])/8;
        end
        function dN = dN_H8(ksi,eta,zeta)
            % dN = dN_H8(ksi,eta,zeta)
            % N [3 x 8][Float]: Shape Functions Derivatives for H8 element
            % Follows anti-clockwise node-numbering convention:
            % (ksi,eta,zeta) = [(-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1)
            %                   (-1,-1,1),  (1,-1,1),  (1,1,1),  (-1,1,1)]
            dN = [];
            for k = [-1 1]
                for j = [-1 1]
                    for i = [-1 1]
                        aux0 = [i*(1+j*eta)*(1+k*zeta) j*(1+i*ksi)*(1+k*zeta) k*(1+j*eta)*(1+i*ksi)]';
                        dN = [dN aux0];
                    end
                end
            end
            dN = dN(:,[1 2 4 3 5 6 8 7])/8;
        end
        function Ndevsparse =  dN_sparse(xi,eta,mu)
            AUX = Element.dN_H8(xi,eta,mu);
            Ndevsparse = [];
            for i = 1:8
                aux0 = AUX(:,i);
                aux1 = [aux0 zeros(3,2)];
                aux2 = [zeros(3,1) aux0 zeros(3,1)];
                aux3 = [zeros(3,2) aux0];
                aux0 = [aux1;aux2;aux3];
                Ndevsparse = [Ndevsparse aux0];
            end
        end
        %% 2D Elements
        %% Q4 Methods
        function N = N_Q4(ksi,eta)
            % N = N_Q4(ksi,eta)
            % N [1 x 4][Float]: Shape Functions for Q4 element
            % Follows anti-clockwise node-numbering convention:
            % (ksi,eta) = [(-1,-1), (1,-1), (1,1), (-1,1)]
            N4 = 0.25*(1 - ksi)*(1 + eta);
            N3 = 0.25*(1 + ksi)*(1 + eta);
            N2 = 0.25*(1 + ksi)*(1 - eta);
            N1 = 0.25*(1 - ksi)*(1 - eta);
            N = [N1 N2 N3 N4];
        end
        function dN = dN_Q4(ksi,eta)
            % dN = dN_Q4(ksi,eta)
            % dN [2 x 4][Float]: Shape functions derivatives
            % dN = [dN_dksi; dN_deta];
            % Follows same numbering convention as N_Q4
            dN = [  % dN_ksi
            -0.25*(1 - eta),    ...
            0.25*(1 - eta),     ...
            0.25*(1 + eta),     ...
            -0.25*(1 + eta)
            % dN_deta
            -0.25*(1 - ksi),    ...
            -0.25*(1 + ksi),    ...
            0.25*(1 + ksi),     ...
            0.25*(1 - ksi) ];
        end
        %% Q8 methods
        function N = N_Q8(ksi,eta)
            % N = N_Q8(ksi,eta)
            % N [1 x 8][Float]: Shape Functions for Q8 element
            % Follows anti-clockwise node-numbering convention:
            % (ksi,eta) = [(-1,-1), (1,-1), (1,1), (-1,1) ...
            %               (0,-1), (1,0),  (0,1), (-1,0)]
            N8 = 0.50*(1 - ksi  )*(1 - eta^2);
            N7 = 0.50*(1 - ksi^2)*(1 + eta  );
            N6 = 0.50*(1 + ksi  )*(1 - eta^2);
            N5 = 0.50*(1 - ksi^2)*(1 - eta  );
            N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8);
            N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7);
            N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6);
            N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8);
            N = [N1 N2 N3 N4 N5 N6 N7 N8];
            
        end
        function dN = dN_Q8(ksi,eta)
            % dN = dN_Q8(ksi,eta)
            % dN [2 x 8][Float]: Shape functions derivatives
            % dN = [dN_dksi; dN_deta];
            % Follows same numbering convention as N_Q8
            dN = [  % dN_ksi
            -0.25*(-1+eta)*(eta+2*ksi),     ...
            -0.25*(-1+eta)*(-eta+2*ksi),    ...
            0.25*(1+eta)*(eta+2*ksi),       ...
            0.25*(1+eta)*(-eta+2*ksi),      ...
            ksi*(-1+eta),                   ...
            -0.5*(-1+eta)*(1+eta),          ...
            -ksi*(1+eta),                   ...
            0.5*(-1+eta)*(1+eta);
            % dN_deta
            -0.25*(-1+ksi)*(ksi+2*eta),     ...
            -0.25*(1+ksi)*(ksi-2*eta),      ...
            0.25*(1+ksi)*(ksi+2*eta),       ...
            0.25*(-1+ksi)*(ksi-2*eta),      ...
            0.5*(-1+ksi)*(1+ksi),           ...
            -(1+ksi)*eta,                   ...
            -0.5*(-1+ksi)*(1+ksi),          ...
            (-1+ksi)*eta ];
        end
        %% Q9 Methods
        function N = N_Q9(ksi,eta)
            % N = N_Q9(ksi,eta)
            % N [1 x 9][Float]: Shape Functions for Q9 element
            % Follows anti-clockwise node-numbering convention:
            % (ksi,eta) = [(-1,-1), (1,-1), (1,1), (-1,1) ...
            %               (0,-1), (1,0),  (0,1), (-1,0), (0,0)]
            N9 =      (1 - ksi^2)*(1 - eta^2);
            N8 = 0.50*(1 - ksi  )*(1 - eta^2) - 0.5*N9;
            N7 = 0.50*(1 - ksi^2)*(1 + eta  ) - 0.5*N9;
            N6 = 0.50*(1 + ksi  )*(1 - eta^2) - 0.5*N9;
            N5 = 0.50*(1 - ksi^2)*(1 - eta  ) - 0.5*N9;
            N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8 + 0.5*N9);
            N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7 + 0.5*N9);
            N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6 + 0.5*N9);
            N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8 + 0.5*N9);
            N  = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
        end
        function dN = dN_Q9(ksi,eta)
            % dN = dN_Q9(ksi,eta)
            % dN [2 x 9][Float]: Shape functions derivatives
            % dN = [dN_dksi; dN_deta];
            % Follows same numbering convention as N_Q9
            dN = [   % dN_ksi
                0.25*eta*(-1+eta)*(2*ksi-1),    ...
                0.25*eta*(-1+eta)*(2*ksi+1),    ...
                0.25*eta*(1+eta)*(2*ksi+1),     ...
                0.25*eta*( 1+eta)*(2*ksi-1),    ...
                -ksi*eta*(-1+eta),              ...
                -1/2*(-1+eta)*(1+eta)*(2*ksi+1),...
                -ksi*eta*(1+eta),               ...
                -1/2*(-1+eta)*(1+eta)*(2*ksi-1),...
                2*ksi*(-1+eta)*(1+eta);
                % dN_eta
                0.25*ksi*(-1+2*eta)*(ksi-1),    ...
                0.25*ksi*(-1+2*eta)*(1+ksi),    ...
                0.25*ksi*(2*eta+1)*(1+ksi),     ...
                0.25*ksi*(2*eta+1)*(ksi-1),     ...
                -0.5*(ksi-1)*(1+ksi)*(-1+2*eta),...
                -ksi*eta*(1+ksi),               ...
                -0.5*(ksi-1)*(1+ksi)*(2*eta+1), ...
                -ksi*eta*(ksi-1),               ...
                2*(ksi-1)*(1+ksi)*eta ];
        end
    end
end