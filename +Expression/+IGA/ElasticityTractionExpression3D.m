classdef ElasticityTractionExpression3D < Expression.IGA.Expression
    %ElasticityBilinearExpression Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        traction_function_
    end
    
    methods
        function this = ElasticityTractionExpression3D()
            this@Expression.IGA.Expression();
        end
        
        function [type, var, basis_id, data] = eval(this, query_unit, differential)
            import Utility.BasicUtility.AssemblyType
            type = AssemblyType.Vector;
            var = {this.test_};
            
            % TODO:In bilinear form consisting of multiplication of derivatives
            % of shape & test, first order derivatives are necessary. Try
            % to improve this process, make it automatic.
            query_unit.query_protocol_{3} = 0;
            
            % Get quadrature
            num_q = query_unit.quadrature_{1};
            qx = query_unit.quadrature_{2};
            qw = query_unit.quadrature_{3};
            
            test_basis = this.test_.basis_data_;
            
            local_matrix = cell(num_q,1);
            
            % loop integration points
            for i = 1 : num_q
                query_unit.query_protocol_{2} = qx(i,:);
                
                % Test query
                test_basis.query(query_unit, []);
                non_zero_id = query_unit.non_zero_id_;
                eval = query_unit.evaluate_basis_;
                             
                % Put non_zero id
                basis_id = {non_zero_id};
                
                % get local mapping
                differential.queryAt(qx(i,:));
                
                [~, J] = differential.jacobian();
                
                x = differential.mapping();
                                                              
                % eval linear form
                N = zeros(3, 3*length(non_zero_id));
                
                id_1 = 1:3:3*length(non_zero_id);
                id_2 = 2:3:3*length(non_zero_id);
                id_3 = 3:3:3*length(non_zero_id);
                
                N(1, id_1) = eval{1};
                
                N(2, id_2) = eval{1};
                
                N(3, id_3) = eval{1};
                
                % add to local matrix
                [t_x, t_y, t_z] = this.traction_function_(x(1), x(2), x(3));
                local_matrix{i} = (N' * [t_x; t_y; t_z]) * qw(i) * J; 
            end
            
            data = local_matrix{1};
            for i = 2 : num_q
                data = data + local_matrix{i};
            end
        end
        
        function setTractionFunction(this, f)
            this.traction_function_ = f;
        end
    end
    
end

