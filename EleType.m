classdef EleType
    enumeration
        Q4, Q8, Q9, AHMAD4, AHMAD8, AHMAD9, H8
    end
    % EleType.Q4, EleType.Q8, EleType.Q9, EleType.H8
    methods (Static)
        function types = Shell_Types()
           types = [EleType.AHMAD4, EleType.AHMAD8, EleType.AHMAD9];
        end
        function bool = is_shell(ele_type)
            shell_types = EleType.Shell_Types();
            bool = false;
            for i = 1:length(shell_types)
                bool = bool || ele_type == shell_types(i);
            end
        end
        
        function N = N_H8(ksi,eta,zeta)
            % N = N_H8(ksi,eta,zeta)
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
    end
end