classdef Operation < handle
    %OPERATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        op_data_
    end
    
    methods
        function this = Operation()
            
        end
        
        function status = setOperator(this, op)
            this.op_data_ = op;
            status = true;
        end
        
        function expression = getExpression(this, method, varargin)
            switch(method)
                case 'FEM'
                    import Expression.FEM.*
                    if strcmp(this.op_data_,'grad_test_dot_grad_var')
                        % Bilinear form
                        if(~isempty(varargin) || ~(length(varargin{1}) < 2))
                            test = varargin{1}{1};
                            variable = varargin{1}{2};
                            expression = BilinearExpression();
                            expression.setTest(test);
                            expression.setVar(variable);
                        else
                            disp('Error <Operation>! - getExpression!');
                            disp('> the bilinear form should input test and variable.');
                        end
                    else
                        expression = Expression();
                    end
                case 'IGA'
                    import Expression.IGA.*
                    switch this.op_data_
                        case 'grad_test_dot_grad_var'
                            % Bilinear form
                            if(~isempty(varargin) || ~(length(varargin{1}) < 2))
                                test = varargin{1}{1};
                                variable = varargin{1}{2};
                                expression = BilinearExpression();
                                expression.setTest(test);
                                expression.setVar(variable);
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the stiffness bilinear form should input test and variable.');
                            end
                        case 'test_dot_var'
                            % Bilinear form for mass matrix
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                variable = varargin{1}{2};
                                penalty_parameter = varargin{1}{3};
                                
                                expression = MassExpression();
                                expression.setTest(test);
                                expression.setVar(variable);
                                if ~isempty(penalty_parameter)
                                    expression.setPenaltyParameter(penalty_parameter);
                                end
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the mass bilinear form should input test and variable.');
                            end
                        case 'test_dot_f'
                            % Linear form for source term
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                source_function = varargin{1}{2};
                                penalty_parameter = varargin{1}{3};
                                
                                expression = LinearExpression();
                                expression.setTest(test);
                                expression.setSourceFunction(source_function);
                                if ~isempty(penalty_parameter)
                                    expression.setPenaltyParameter(penalty_parameter);
                                end
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the linear form should input test and variable.');
                            end
                        case 'delta_epsilon_dot_sigma'
                            % Bilinear form for linear elasticity
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                variable = varargin{1}{2};
                                constitutive_law = varargin{1}{3};
                                
                                expression = ElasticityBilinearExpression(constitutive_law);
                                expression.setTest(test);
                                expression.setVar(variable);
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the stiffness bilinear form should input test and variable.');
                            end
                        case 'ElasticityBilinearExpression3D'
                            % Bilinear form for linear elasticity
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                variable = varargin{1}{2};
                                constitutive_law = varargin{1}{3};
                                
                                expression = ElasticityBilinearExpression3D(constitutive_law);
                                expression.setTest(test);
                                expression.setVar(variable);
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the stiffness bilinear form should input test and variable.');
                            end    
                        case 'delta_Tr(epsilon)_T'
                            % Linear form for linear elasticity
                            if(~isempty(varargin) || ~(length(varargin{1}) < 2))
                                test = varargin{1}{1};
                                thermal_function = varargin{1}{2};
                                
                                expression = ElasticityLinearExpression();
                                expression.setTest(test);
                                expression.setSourceFunction(thermal_function);
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the stiffness bilinear form should input test and variable.');
                            end   
                        case 'laplace_nitsche_dirichlet_lhs_term'
                            % Laplace Dirichlet boundary term using Nitsche's method
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                variable = varargin{1}{2};
                                penalty_parameter = varargin{1}{3};
                                
                                expression = NitscheLhsExpression();
                                expression.setTest(test);
                                expression.setVar(variable);
                                
                                if ~isempty(penalty_parameter)
                                    expression.setPenaltyParameter(penalty_parameter);
                                end
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the Nitsche lhs term should input test, source function and beta.');
                            end
                        case 'laplace_nitsche_dirichlet_rhs_term'
                            % Laplace Dirichlet boundary term using Nitsche's method
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                source_function = varargin{1}{2};
                                penalty_parameter = varargin{1}{3};
                                
                                expression = NitscheRhsExpression();
                                expression.setTest(test);                               
                                expression.setSourceFunction(source_function);
                                
                                if ~isempty(penalty_parameter)
                                    expression.setPenaltyParameter(penalty_parameter);
                                end
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the Nitsche lhs term should input test, variable and beta.');
                            end
                        case 'laplace_nitsche_interface_lhs_term'
                            % Laplace Interface term using Nitsche's method
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test_1 = varargin{1}{1};
                                variable_1 = varargin{1}{2};
                                test_2 = varargin{1}{3};
                                variable_2 = varargin{1}{4};
                                penalty_parameter = varargin{1}{5};
                                
                                expression = NitscheInterfaceExpression();
                                expression.setTest({test_1, test_2});
                                expression.setVar({variable_1, variable_2});
                                
                                if ~isempty(penalty_parameter)
                                    expression.setPenaltyParameter(penalty_parameter);
                                end
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the Nitsche interface term should input two test, two variable functions and beta.');
                            end 
                        case '3D_nonlinear_elasticity'    
                            % Tangent matrix and resiual fores for nonlinear elasticity
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                variable = varargin{1}{2};
                                constitutive_law = varargin{1}{3};
                                
                                expression = ElasticityNonlinearExpression3D(constitutive_law);
                                expression.setTest(test);
                                expression.setVar(variable);
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the stiffness bilinear form should input test and variable.');
                            end
                        case '2D_nonlinear_elasticity'    
                            % Tangent matrix and resiual fores for nonlinear elasticity
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                variable = varargin{1}{2};
                                constitutive_law = varargin{1}{3};
                                
                                expression = ElasticityNonlinearExpression2D(constitutive_law);
                                expression.setTest(test);
                                expression.setVar(variable);
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the stiffness bilinear form should input test and variable.');
                            end
                        case 'elasticity_nitsche_dirichlet_lhs_term'
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                variable = varargin{1}{2};
                                constitutive_law = varargin{1}{3};
                                penalty_parameter = varargin{1}{4};
                                                              
                                expression = ElasticityNitscheLhsExpression(constitutive_law);
                                expression.setTest(test);
                                expression.setVar(variable);
                                
                                if ~isempty(penalty_parameter)
                                    expression.setPenaltyParameter(penalty_parameter);
                                end
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the Nitsche lhs term should input test, source function and beta.');
                            end    
                        case 'elasticity_nitsche_dirichlet_rhs_term'
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                source_function = varargin{1}{2};
                                penalty_parameter = varargin{1}{3};
                                
                                expression = ElasticityNitscheRhsExpression();
                                expression.setTest(test);                               
                                expression.setSourceFunction(source_function);
                                
                                if ~isempty(penalty_parameter)
                                    expression.setPenaltyParameter(penalty_parameter);
                                end
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the Nitsche rhs term should input test, source and beta.');
                            end
                        case 'elasticity_applied_traction'
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                traction_function = varargin{1}{2};
                                
                                expression = ElasticityTractionExpression();
                                expression.setTest(test); 
                                expression.setTractionFunction(traction_function);
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the Nitsche lhs term should input test and traction.');
                            end  
                        case 'ElasticityTractionExpression3D'
                            if(~isempty(varargin) || ~(length(varargin{1}) < 3))
                                test = varargin{1}{1};
                                traction_function = varargin{1}{2};
                                
                                expression = ElasticityTractionExpression3D();
                                expression.setTest(test); 
                                expression.setTractionFunction(traction_function);
                            else
                                disp('Error <Operation>! - getExpression!');
                                disp('> the Traction term should input test and traction.');
                            end    
                        otherwise
                            expression = Expression();
                    end   
                otherwise
                    expression = [];
                    disp('Error <Operation>! - getExpression!');
                    disp('> Your input method not supported. ');
            end
        end
    end
    
end

