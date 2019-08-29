classdef ElasticityBilinearExpression3D < Expression.IGA.Expression
    %ElasticityBilinearExpression Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        constitutive_law_
    end
    
    methods
        function this = ElasticityBilinearExpression3D(constitutive_law)
            this@Expression.IGA.Expression();
            this.constitutive_law_ = constitutive_law;
        end
        
        function [type, var, basis_id, data] = eval(this, query_unit, differential)
            import Utility.BasicUtility.AssemblyType
            type = AssemblyType.Matrix;
            var = {this.test_; this.var_};
            
            % TODO:In bilinear form consisting of multiplication of derivatives
            % of shape & test, first order derivatives are necessary. Try
            % to improve this process, make it automatic.
            query_unit.query_protocol_{3} = 1;
            
            % Get quadrature
            num_q = query_unit.quadrature_{1};
            qx = query_unit.quadrature_{2};
            qw = query_unit.quadrature_{3};
            
            test_basis = this.test_.basis_data_;
            
            local_matrix = cell(num_q,1);
            
            D_matrix = this.constitutive_law_.materialMatrix();
            % loop integration points
            for i = 1 : num_q
                query_unit.query_protocol_{2} = qx(i,:);
                
                % Test query
                test_basis.query(query_unit, []);
                non_zero_id = query_unit.non_zero_id_;
                eval = query_unit.evaluate_basis_;
              
                % Put non_zero id
                basis_id = {non_zero_id, non_zero_id};
                
                % get local mapping
                differential.queryAt(qx(i,:));

                [dx_dxi, J] = differential.jacobian();
                                              
                % eval basis derivative with x
                d_N_dx = dx_dxi \ eval{2};
                
                % eval bilinear form
                B_matrix = zeros(6, 3*length(non_zero_id));
                
                id_1 = 1:3:3*length(non_zero_id);
                id_2 = 2:3:3*length(non_zero_id);
                id_3 = 3:3:3*length(non_zero_id);
                
                B_matrix(1, id_1) = d_N_dx(1,:);
                
                B_matrix(2, id_2) = d_N_dx(2,:);
                
                B_matrix(3, id_3) = d_N_dx(3,:);
                
                B_matrix(4, id_1) = d_N_dx(2,:);
                B_matrix(4, id_2) = d_N_dx(1,:);
                      
                B_matrix(5, id_2) = d_N_dx(3,:);
                B_matrix(5, id_3) = d_N_dx(2,:);
                
                B_matrix(6, id_3) = d_N_dx(1,:);
                B_matrix(6, id_1) = d_N_dx(3,:);
                
                % add to local matrix
                local_matrix{i} = (B_matrix' * D_matrix * B_matrix).* qw(i) * J; 
            end
            
            data = local_matrix{1};
            for i = 2 : num_q
                data = data + local_matrix{i};
            end
        end

    end
    
end

