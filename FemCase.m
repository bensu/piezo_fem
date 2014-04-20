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
    end
end