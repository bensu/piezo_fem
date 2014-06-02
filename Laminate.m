classdef Laminate
    % class Laminate
    % Data structure holding a list of materials and their positions
    properties
        mat_list
        t_list
    end
    methods 
        function obj = Laminate(mat_list,t_list)
            require(all(isa(mat_list,'Material')), ...
                'ArgumentError: mat_list should be a list of materials');
            require(all(isnumeric(t_list)), ...
                'ArgumentError: t_list should be numeric');
            require(all(size(mat_list)==size(t_list)), ...
                'ArgumentError: mat_list and t_list should have the same size');
            obj.mat_list = mat_list;
            obj.t_list   = t_list;
        end
    end
    methods
        function mat_out = material(lam,zeta)
            % mat_out = material(lam,zeta)
            % mat_out [Material]: Corresponds to layer in location zeta
            % zeta [Float]: between -1 and 1, local coordinate
            % Maybe it should be a thickness value!
            zeta_vals = lam.zeta_values;
            mat_out = lam.mat_list(find(zeta <= zeta_vals,1,'last'));            
        end
        function zeta_vals = zeta_values(lam)
            % zeta_vals = zeta_values(lam)
            % Table to interpolate from zeta values to thickness values
            t_sum = sum(lam.t_list);    % Top surface location
            zeta_vals = interp1([0 t_sum/2 t_sum],[-1 0 1], ...
                                    cumsum(lam.t_list),'linear','extrap');
        end
    end
end