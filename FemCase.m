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
            obj.bc = CompoundFunction(true,mesh_in.n_nodes, ...
            	physics_in.dofs_per_node,mesh_in.n_ele,physics_in.dofs_per_ele);
            obj.loads = CompoundFunction(0,mesh_in.n_nodes, ...
            	physics_in.dofs_per_node,mesh_in.n_ele,physics_in.dofs_per_ele);
            obj.dis = CompoundFunction(0,mesh_in.n_nodes, ...
            	physics_in.dofs_per_node,mesh_in.n_ele,physics_in.dofs_per_ele);
            obj.reactions = CompoundFunction(0,mesh_in.n_nodes, ...
            	physics_in.dofs_per_node,mesh_in.n_ele,physics_in.dofs_per_ele);
        end
        function solve(fem)
            % Side Effects
            % Solves the problem set by physics.k and the mesh.
            % BC and Loads to DOF form
            L = fem.loads.all_dofs;
            F = fem.bc.all_dofs;
            % Create Stiffness
            S = fem.mesh.assembly(fem.physics.dofs_per_node, ...
                                  fem.physics.dofs_per_ele,  ...
                                  fem.physics.k);
            D = zeros(size(S,1),1);
            D(F) = S(F,F) \ L(F);
            fem.dis.dof_list_in(D);
            R = zeros(size(S,1),1);
%             R(F) = S(F,F) * D(F);
            R = S * D;
            fem.reactions.dof_list_in(R);
            edge = fem.mesh.find_nodes(@(x,y,z) (abs(x)<1e-5));
            fem.reactions.node_vals.vals(edge,3)
            sum(abs(fem.reactions.node_vals.vals))
            sum(fem.loads.node_vals.vals)
        end
    end
end