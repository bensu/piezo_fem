classdef CompoundFunction < hgsetget & handle
    properties
        ele_vals
        node_vals
    end
    properties (Dependent)
        n_ele
        n_nodes
        dofs_per_ele
        dofs_per_node
        n_ele_dofs
        n_node_dofs
        all_dofs
    end
    methods
        function obj = CompoundFunction(filler,total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in)
            node_vals_in = VectorFunction(filler,total_number_of_nodes_in,dofs_per_node_in);
            ele_vals_in =  VectorFunction(filler,total_number_of_elements_in,dofs_per_element_in);
            obj.set('ele_vals',ele_vals_in);
            obj.set('node_vals',node_vals_in);
        end
        
        %% I/O
        function dofs = get.all_dofs(cf)
            dofs = [cf.get('node_vals').all_dofs;
                       cf.get('ele_vals').all_dofs];
            
        end
        function dof_list_in(compound,dofs)
            nodes_dofs = dofs(1:compound.n_node_dofs);
            element_dofs = dofs((compound.n_node_dofs+1):end);
            assert(length(element_dofs)==compound.n_ele_dofs);
            compound.get('node_vals').dof_list_in(nodes_dofs);
            compound.get('ele_vals').dof_list_in(element_dofs);
        end        
        function fun_out = node_function(compound)
            fun_out = compound.get('node_vals').get('component_function');
        end
        function fun_out = element_function(compound)
            fun_out = compound.get('ele_vals').get('component_function');
        end
        
        function clear(compound)
            compound.get('node_vals').clear
            compound.get('ele_vals').clear
        end
        function dofs = non_zero_index(compound)
            % dofs = non_zero_index(fun)
            % Returns the logical index of all the values that are not zero
            dofs = [compound.node_vals.non_zero_index
                    compound.ele_vals.non_zero_index];
        end
        
        %% Helpers
        function num = get.n_nodes(compound)
            % Number of nodes
            num = compound.get('node_vals').get('total_number_of_components');            
        end
        function num = get.n_ele(compound)
            % Number of Elements
            num = compound.get('ele_vals').get('total_number_of_components');            
        end
        function num = get.dofs_per_node(compound)
            % Dofs per node
            num = compound.get('node_vals').get('dofs_per_component');            
        end
        function num = get.dofs_per_ele(compound)
            % Dofs per element
            num = compound.get('node_element').get('dofs_per_component');           
        end
        function num = get.n_node_dofs(compound)
            % Total number of node dofs
            num = compound.get('node_vals').number_of_dofs;            
        end
        function num = get.n_ele_dofs(compound)
            % Total number of ele dofs
            num = compound.get('ele_vals').number_of_dofs;            
        end
    end                    
end