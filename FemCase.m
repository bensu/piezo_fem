classdef FemCase < handle
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
            obj.bc = CompoundFunction(false,mesh_in.n_nodes, ...
            	physics_in.dofs_per_node,mesh_in.n_ele,physics_in.dofs_per_ele);
            obj.loads = CompoundFunction(0,mesh_in.n_nodes, ...
            	physics_in.dofs_per_node,mesh_in.n_ele,physics_in.dofs_per_ele);
            obj.dis = CompoundFunction(0,mesh_in.n_nodes, ...
            	physics_in.dofs_per_node,mesh_in.n_ele,physics_in.dofs_per_ele);
            obj.reactions = CompoundFunction(0,mesh_in.n_nodes, ...
            	physics_in.dofs_per_node,mesh_in.n_ele,physics_in.dofs_per_ele);
        end
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
        function [V,D] = eigen_values(fem)
        % Solves the eigen value problem for K and M.
            require(~isempty(fem.physics.m), ...
                'Error: physics needs a mass matrix');
            F = ~fem.bc.all_dofs;
            S = fem.S;
            M = fem.M;
            [V,D] = eigs(M(F,F),S(F,F));
            
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
%                 cond(S(F,F))
%                 cond(S_scaled(F,F))
%                 size(S_scaled(F,F))
%                 rank(S_scaled(F,F))
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
    end
end