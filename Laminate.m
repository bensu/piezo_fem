classdef Laminate
    % class Laminate
    % Data structure holding a list of materials and their positions
    % Represents layers by 
    properties
        mat_list
        t_list
    end
    properties (Dependent)
        zeta_t      % Value of thickness distance in zeta coordinates, positive
        t_m         % 
        zeta_m
        n_layers
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
        function C = PiezoMatrix(lam,zeta)
            % C = PiezoMatrix(lam,zeta)
            % Get the constitutive matrix for a certain layer
            % Get Layer's contribution
            n_l = lam.n_layers;
            layer = lam.material(zeta);
            l = lam.mat_num(zeta);
            elastic = Physics.ElasticShell(layer);
            piezo = -layer.D(:,[1 2 4 5 6]);
            electric = -diag(layer.e);
            % Put the layer's electric contribution in the general matrix
            all_piezo = zeros(n_l*size(piezo,1),size(piezo,2));
            all_piezo(index_range(size(piezo,1),l),:) = piezo;
            all_electric = zeros(size(electric,1)*n_l);
            e_index = index_range(size(electric,2),l);
            all_electric(e_index,e_index) = electric;
            C = [   elastic     all_piezo';
                all_piezo   all_electric];
        end
        function mat_out = material(lam,zeta)
            % mat_out = material(lam,zeta)
            % mat_out [Material]: Corresponds to layer in location zeta
            % zeta [Float]: between -1 and 1, local coordinate
            % Maybe it should be a thickness value!
            mat_out = lam.mat_list(lam.mat_num(zeta));            
        end
        function l = mat_num(lam,zeta)
            % l = mat_num(lam,zeta)
            % Finds the layer number that corresponds to zeta
            zeta_vals = lam.zeta_top;
            l = find(zeta <= zeta_vals,1,'first');
        end
        function [zeta_p, zeta_w] = quadrature(laminate,order)
            % [zeta_p, zeta_w] = quadrature(laminate,order)
            % Generates all the gauss points and weights for a piecewise
            % integration along the laminate.
            zeta_top = [-1 laminate.zeta_top];
            zeta_p = [];
            zeta_w = [];
            for i = 2:length(zeta_top)
                [g_p,g_w] = Integral.lgwt(order,zeta_top(i-1),zeta_top(i));
                zeta_p    = [zeta_p g_p];
                zeta_w    = [zeta_w g_w];
            end
        end
        function zeta_vals = zeta_top(lam)
            % zeta_vals = zeta_top(lam)
            % Transforms the thickness coordinate of each layer's top surface to
            % zeta coordinate
            zeta_vals = lam.t_to_z(cumsum(lam.t_list));
        end
        function z = t_to_z(lam,t)
            % Transforms a value in thickness coordinate to its equivalent
            % zeta local coordinate
            t_sum = sum(lam.t_list);    % Top surface location
            z = interp1([0 t_sum/2 t_sum],[-1 0 1],t,'linear','extrap');
        end
        function t = z_to_t(lam,z)
            % t = z_to_t(lam,z)
            % Inverse of t_to_z
            t_sum = sum(lam.t_list);    % Top surface location
            t = interp1([-1 0 1],[0 t_sum/2 t_sum],z,'linear','extrap');
        end
        %% Dependent Properties
        function zeta_t = get.zeta_t(lam)
            zeta_t = lam.t_to_z(lam.t_list) + 1;
        end
        function zeta_m = get.zeta_m(lam)
            % zeta_m = mid_surface_locations(lam)
            % mid_surface_locations
            % zeta_m [n_layers x 1][Float]: Returns the location of each
            % layer's midsurface in local coordinates
            zeta_m = lam.t_to_z(lam.t_m);
        end
        function t_m = get.t_m(lam)
            % t_m = mid_surface_locations(lam)
            % mid_surface_locations
            % t_m [n_layers x 1][Float]: Returns the location of each
            % layer's midsurface.
            t_m = lam.t_list/2 + cumsum(lam.t_list) - lam.t_list(1);
        end
        function n = get.n_layers(lam)
            % n = get.n_layers(lam)
            % Number of layers
            n = length(lam.t_list);
        end
    end
end