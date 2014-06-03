classdef VectorFunction < hgsetget
    properties
        vals
        dofs_per_val
        n_vals
        zero_val
    end
    properties (Dependent)
        all_dofs
    end
    methods
        function obj = VectorFunction(filler,n_vals_in, dofs_per_val_in)
            % Always stored in val and element function
            val_list_in = repmat(filler,n_vals_in,dofs_per_val_in);
            obj.set('dofs_per_val',dofs_per_val_in)
            obj.set('vals',val_list_in)
            obj.set('n_vals',n_vals_in)
            obj.set('zero_val',filler);
        end
        
        %% I/O
        function dof_list_in(fun,component_list)
            fun_in = reshape(component_list,fun.get('dofs_per_val'), ...
                        fun.get('n_vals'))';
            fun.set('vals',fun_in);
        end
        function vectors_out = component_list(fun)
            vectors_out = fun.get('vals');
        end
        function dofs_out = get.all_dofs(fun)
            dofs_out = reshape(fun.get('vals')',[],1);
        end
        function dofs = non_zero_index(fun)
            % dofs = non_zero_index(fun)
            % Returns the logical index of all the values that are not zero
            dofs = (fun.all_dofs ~= fun.zero_val);
        end
        
        %% Edit
        function set_val(vf,places,value)
            % Sets all the values where places is true to value
            for v = 1:vf.n_vals
                if places(v)
                    vf.vals(v,:) = value;
                end
            end
        end
        function clear(vector_fun)
            filler = vector_fun.get('zero_val');
            vector_fun.set('vals',repmat(filler, ...
                vector_fun.get('n_vals'), ...
                vector_fun.get('dofs_per_val')));
        end
        function edit_component_by_id(vector,component_id,vector_in)
            require(length(vector_in)==vector.get('dofs_per_val'), ...
                    'ArgumentError: Wrong vector_in dimension')
            fun_list = vector.get('vals');
            fun_list(component_id,:) = vector_in;
            vector.set('val',fun_list);
        end
        function edit_component_part_by_id(vector,component_id,part,scalar_in)
            % Only for scalar values
            require(part<=vector.get('dofs_per_val'),  ...
                'ArgumentError: Component too big')
            fun_list = vector.get('vals');
            fun_list(component_id,part) = scalar_in;
            vector.set('vals',fun_list);
        end     
        
        %% Helpers
        function num = number_of_dofs(vector_fun)
            num = vector_fun.get('dofs_per_val')* ...
                        vector_fun.get('n_vals');
            if isempty(num)
                num = 0;
            end
        end
    end
end