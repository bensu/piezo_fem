classdef VectorFunction < hgsetget
    properties
        component_function
        dofs_per_component
        total_number_of_components
        zero_val
    end            
    methods
        function obj = VectorFunction(filler, ... 
                        total_number_of_components_in, dofs_per_component_in)
            % Always stored in component_function and element function
            component_list_in = repmat(filler,total_number_of_components_in, ...
                            dofs_per_component_in);
            obj.set('dofs_per_component',dofs_per_component_in)
            obj.set('component_function',component_list_in)
            obj.set('total_number_of_components',total_number_of_components_in)
            obj.set('zero_val',filler);
        end
        
        %% I/O
        function dof_list_in(fun,component_list)
            fun_in = reshape(component_list,fun.get('dofs_per_component'), ...
                        fun.get('total_number_of_components'))';
            fun.set('component_function',fun_in);
        end
        function vectors_out = component_list(fun)
            vectors_out = fun.get('component_function');
        end
        function dofs_out = all_dofs(fun)
            dofs_out = reshape(fun.get('component_function')',[],1);
        end
        
        %% Edit
        function clear(vector_fun)
            filler = vector_fun.get('zero_val');
            vector_fun.set('component_function',repmat(filler, ...
                vector_fun.get('total_number_of_components'), ...
                vector_fun.get('dofs_per_component')));
        end
        
        function edit_component_by_id(vector,component_id,vector_in)
            require(length(vector_in)==vector.get('dofs_per_component'), ...
                    'Wrong vector_in dimension')
            fun_list = vector.get('component_function');
            fun_list(component_id,:) = vector_in;
            vector.set('component_function',fun_list);
        end
        function edit_component_part_by_id(vector,component_id,part,scalar_in)
            % Only for scalar values
            require(part<=vector.get('dofs_per_component'), 'Wrong component too big')
            fun_list = vector.get('component_function');
            fun_list(component_id,part) = scalar_in;
            vector.set('component_function',fun_list);
        end     
        
        %% Helpers
        function num = number_of_dofs(vector_fun)
            num = vector_fun.get('dofs_per_component')* ...
                        vector_fun.get('total_number_of_components');
            if isempty(num)
                num = 0;
            end
        end
    end
end