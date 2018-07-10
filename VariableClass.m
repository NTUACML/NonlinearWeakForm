classdef VariableClass
    %VARIABLECLASS Summary of this class goes here
    %   Detailed explanation goes here
    properties
        data_           % variable data
        name_           % variable name
        dim_            % variable dimension
        data_number_    % total variable data number
    end
    
    methods
        % constructor
        function this = VariableClass(name, dim, basis_num)
            this.name_ = name;
            this.dim_ = dim;
            if nargin < 3
                this.data_ =  [];
                this.data_number_ = 0;
            else
                this.data_ =  zeros(dim * basis_num, 1);
                this.data_number_ = dim * basis_num;
            end
        end
        % function to obtaine data component
        function data_out = data_component(this, component_dim)
            data_out = this.data_(component_dim:this.dim_:this.data_number_);
        end
    end
end
