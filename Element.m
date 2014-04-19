classdef Element
    properties
        coords
        normals
    end
    methods
        function obj = Element(coords,normals)
            require(size(coords,1)==4, ...
                'ArgumentError: only 4 nodes');
            require(size(coords)==size(normals), ...
                'ArgumentError: coords and normals should have same size');
            obj.coords = coords;
            obj.normals = normals;
        end
        function jac_out = jacobian(element,xi,eta,mu)
            % Computes the jacobian for the ShellQ4 element
            dNdxi = Element.dNdxi_Q4(xi,eta);
            dNdeta = Element.dNeta_Q4(xi,eta);
            N = Element.N_Q4(xi,eta);
            % This ASSUMES element.normals has a certain shape
            jac_out = [dNdxi;dNdeta;zeros(1,4)]*element.coords ...
                + [mu*dNdxi;mu*dNdeta;N]*element.normals/2;
        end
    end
    methods (Static)
        function N_out = N_ShellQ4(element,xi,eta,mu)
            % Not really used!!!
            require(isnumeric(mu), ...
                'ArgumentError: xi, eta, and mu should be numeric')
            require(-1<=mu && mu<=1, ... 
                'ArgumentError: mu should is not -1<=mu<=1')
            % Need to check the way it works with Jacobian
            N_out = Element.N_Q4(xi,eta)*mu*element.normals;
        end
        function N_out = N_Q4(xi,eta)
            %  Notes:
            %     1st node at (-1,-1), 3rd node at (-1,1) 
            %     4th node at (1,1), 2nd node at (1,-1)
            require(isnumeric([xi eta]), ...
                'ArgumentError: Both xi and eta should be numeric')
            require(-1<=xi && xi<=1, ...
                'ArgumetnError: xi should be -1<=xi<=1')
            require(-1<=eta && eta<=1, ...
                'ArgumetnError: eta should be -1<=eta<=1')
            N_out = zeros(1,4);
            N_out(1) = 0.25*(1-xi)*(1-eta);
            N_out(3) = 0.25*(1+xi)*(1-eta);
            N_out(2) = 0.25*(1-xi)*(1+eta);
            N_out(4) = 0.25*(1+xi)*(1+eta);
        end
        function dNdxiQ4_out = dNdxi_Q4(xi,eta)
            % derivatives
            require(isnumeric([xi eta]), ...
                'ArgumentError: Both xi and eta should be numeric')
            require(-1<=xi && xi<=1, ...
                'ArgumetnError: xi should be -1<=xi<=1')
            require(-1<=eta && eta<=1, ...
                'ArgumetnError: eta should be -1<=eta<=1')
            dNdxiQ4_out = zeros(1,4);
            dNdxiQ4_out(1) = -0.25*(1-eta);
            dNdxiQ4_out(2) = 0.25*(1-eta);
            dNdxiQ4_out(3) = -0.25*(1+eta);
            dNdxiQ4_out(4) = 0.25*(1+eta);
        end
        function dNdetaQ4_out = dNeta_Q4(xi,eta)
            % derivatives
            require(isnumeric([xi eta]), ...
                'ArgumentError: Both xi and eta should be numeric')
            require(-1<=xi && xi<=1, ...
                'ArgumetnError: xi should be -1<=xi<=1')
            require(-1<=eta && eta<=1, ...
                'ArgumetnError: eta should be -1<=eta<=1')
            dNdetaQ4_out = zeros(1,4);
            dNdetaQ4_out(1) = -0.25*(1-xi);
            dNdetaQ4_out(2) = -0.25*(1+xi);
            dNdetaQ4_out(3) = 0.25*(1-xi);
            dNdetaQ4_out(4) = 0.25*(1+xi);
        end
    end
end
            