classdef FemCase < handle
    % FemCase Class
    % Works a the collection of data structures that represent a FEM
    % problem. It includes a mesh, the physical problem, it's boundary
    % conditions and the results.
    % Main Methods:
    %   - solve
    %   - eigen_values
    properties
        mesh
        physics
        loads
        bc
        dis
        reactions
    end
    methods
        function obj = FemCase(mesh_in,physics_in)
            require(isa(mesh_in,'Mesh'), ...
                'ArgumentError: mesh_in not type Mesh');
            require(isa(physics_in,'Physics'), ...
                'ArgumentError: mesh_in not type Physics');
            obj.mesh = mesh_in;
            obj.physics = physics_in;
            % Create the rest of the Properties
            
            obj.bc = obj.compound_function(false);
            obj.loads = obj.compound_function(0);
            obj.dis = obj.compound_function(0);
            obj.reactions = obj.compound_function(0);
        end
        function [V,D] = eigen_values(fem,mode_number)
            % [V,D] = eigen_values(fem,mode_number)
            % Solves the eigen value problem for K and M.
            % V [n_dof x n_dof][Float]: Mode Shape Matrix
            % D [n_dof x n_dof][Float]: Diagonal matrix with omega^2
            % mode_number [Int]: Number of modes to be calculated
            require(~isempty(fem.physics.m), ...
                'Error: physics needs a mass matrix');
            F = ~fem.bc.all_dofs;
            S = fem.S;
            M = fem.M;
            [V1,D] = eigs(S(F,F),M(F,F),mode_number,'SM');
            V = fem.compound_function(0);
            aux = V.all_dofs;
            aux(F) = V1(:,3);
            V.dof_list_in(aux);
        end
        function solve(fem)
            % Side Effects
            % Solves the problem set by physics.k and the mesh.
            % BC and Loads to DOF form
            L = fem.loads.all_dofs;
            F = ~fem.bc.all_dofs;
            % Create Stiffness
            S = fem.mesh.assembly_matrix(fem.physics.dofs_per_node, ...
                                         fem.physics.dofs_per_ele,  ...
                                         fem.physics.k);
            D = zeros(size(S,1),1);
            
            % For unit reasons (i.e. piezoelectric constant, some terms in 
            % the stiffness matrix (S) are several orders of magnitude 
            % larger than others. It is desirable to improve cond(S) before
            % solving the system.
            if (condest(S(F,F)) > 1e8)
                P = diag(sqrt(diag(S)));
                inv_P = inv(P);
                S_scaled = inv_P*S*inv_P;
                L_scaled = inv_P*L;
                D_scaled = S_scaled(F,F) \ L_scaled(F);
                D(F) = inv_P(F,F)*D_scaled;
                % Change to true if output is wanted
                if (false)
                    cond(S(F,F))
                    cond(S_scaled(F,F))
                    size(S_scaled(F,F))
                    rank(S_scaled(F,F))
                end
            else
                D(F) = S(F,F) \ L(F);
            end

            fem.dis.dof_list_in(D);
            fem.reactions.dof_list_in(S*D);
            edge = fem.mesh.find_nodes(@(x,y,z) (abs(x)<1e-5));
            fem.reactions.node_vals.vals(edge,3);
            sum(abs(fem.reactions.node_vals.vals));
            sum(fem.loads.node_vals.vals);
        end
        %% Wrappers
        function S = S(fem)
        % S = S(fem)
        % Wrapper for Stiffness Matrix
            S = fem.mesh.assembly_matrix(fem.physics.dofs_per_node, ...
                fem.physics.dofs_per_ele,  ...
                fem.physics.k);
        end
        function M = M(fem)
        % M = M(fem)
        % Wrapper for Mass Matrix
            M = fem.mesh.assembly_matrix(fem.physics.dofs_per_node, ...
                fem.physics.dofs_per_ele,  ...
                fem.physics.m);
        end
        function obj = compound_function(fem,filler)
            % Wrapper for class CompoundFunction
            obj = CompoundFunction(filler,fem.mesh.n_nodes, ...
            	fem.physics.dofs_per_node,fem.mesh.n_ele,fem.physics.dofs_per_ele);
        end
        function plot(fem)
            if ~isempty(fem.dis) && ~isempty(fem.reactions)
                fem.mesh.plot_displacements(100,fem.dis.node_vals.vals(:,1:3))
                hold on
                fem.mesh.plot_vector_field( ...
                    fem.reactions.node_vals.vals(:,1:3),'k');
            else
                hold on
                fem.mesh.plot
            end
            fem.mesh.plot_vector_field(fem.bc.node_vals.vals(:,1:3),'b')
            fem.mesh.plot_vector_field(fem.loads.node_vals.vals(:,1:3),'r')
            hold off
        end
    end
end