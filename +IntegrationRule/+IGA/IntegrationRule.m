classdef IntegrationRule < IntegrationRule.IntegrationRuleBase
    %INTEGRATIONRULE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function this = IntegrationRule(patch, expression)
            this@IntegrationRule.IntegrationRuleBase(patch, expression);
        end
        
        function status = generate(this, generate_parameter)
            import Utility.BasicUtility.Region
            if(this.integral_patch_.region_ == Region.Domain)
                status = this.generateDomainIntUnit(generate_parameter);
            elseif(this.integral_patch_.region_ == Region.Boundary)
                status = this.generateBoundaryIntUnit(generate_parameter);
            else
                disp('Error <IntegrationRule>! - generate!');
                disp('> Your integral region is wrong! please check it~');
                status = false;
            end
        end
    end
    
    methods(Access = private)
        function status = generateDomainIntUnit(this, generate_parameter)
            nurbs_data = this.integral_patch_.nurbs_data_;  
            if isempty(generate_parameter)
                generation_method = 'Default';
                number_quad_pnt = ceil((nurbs_data.order_ + 1)*0.5);  
            else
                generation_method = generate_parameter{1};
                number_quad_pnt = generate_parameter{2}*ones(size(nurbs_data.order_));
            end
            
            % IntUnit number
            switch generation_method
                case 'Default'
                    this.generateIntegralUnitByKnotMesh(nurbs_data);
                case 'QuadTree'
                    disp('To be finihed...');
                    %nurbs_mesh = this.generateIntegralUnitByQuadTree(some_stratergy);
            end
            
            
            % generate quadrature rule
            import IntegrationRule.IGA.GaussQuadrature
            for int_unit_id = 1 : this.num_integral_unit_
                % add guass quadrature data
                this.integral_unit_{int_unit_id}.quadrature_ = ...
                    GaussQuadrature.MappingNurbsType2GaussQuadrature(...
                    nurbs_data.type_, this.integral_unit_{int_unit_id}.unit_span_, number_quad_pnt);
                
            end
            status = true;
        end
        
        function status = generateBoundaryIntUnit(this, generate_parameter)
            if isa(this.integral_patch_, 'Utility.BasicUtility.InterfacePatch')            
                nurbs_data = this.integral_patch_.master_patch_.nurbs_data_;  
            else
                nurbs_data = this.integral_patch_.nurbs_data_;
            end
            
            if isempty(generate_parameter)
                generation_method = 'Default';
                number_quad_pnt = ceil((nurbs_data.order_ + 1)*0.5);  
            else
                generation_method = generate_parameter{1};
                number_quad_pnt = generate_parameter{2};
            end
                        
            % IntUnit number
            switch generation_method
                case 'Default'
                    this.generateIntegralUnitByKnotMesh(nurbs_data);
                case 'QuadTree'
                    disp('To be finihed...');
                    %nurbs_mesh = this.generateIntegralUnitByQuadTree(some_stratergy);
            end
            
            
            % generate quadrature rule
            import IntegrationRule.IGA.GaussQuadrature
            for int_unit_id = 1 : this.num_integral_unit_
                % add guass quadrature data
                integral_unit = this.integral_unit_{int_unit_id};
                integral_unit.quadrature_ = ...
                    GaussQuadrature.MappingNurbsType2GaussQuadrature(...
                    nurbs_data.type_, this.integral_unit_{int_unit_id}.unit_span_, number_quad_pnt);                
            end
            
            % if the patch is a interface patch, the integration units are
            % transformed into those consisted of only one quadrature
            % point. However, this is not a good way to do that. The
            % generating procesure needs to be re-designed.
            if isa(this.integral_patch_, 'Utility.BasicUtility.InterfacePatch')            
                total_num_quad_point = this.num_integral_unit_*number_quad_pnt;
                qx = zeros(total_num_quad_point, this.integral_patch_.dim_);
                wx = zeros(total_num_quad_point, 1);
                
                for i = 1:this.num_integral_unit_
                    for j = 1:number_quad_pnt
                        qx((i-1)*number_quad_pnt+j,:) = this.integral_unit_{i}.quadrature_{2}(j,:);
                        wx((i-1)*number_quad_pnt+j) = this.integral_unit_{i}.quadrature_{3}(j);
                    end
                end
                
                new_integral_unit = cell(total_num_quad_point,1);
                
                import IntegrationRule.IGA.IntegralUnit
                for i = 1:this.num_integral_unit_
                    for j = 1:number_quad_pnt
                        cnt = (i-1)*number_quad_pnt+j;
                        new_integral_unit{cnt} = IntegralUnit(this.integral_patch_, this.integral_unit_{i}.unit_span_);
                        new_integral_unit{cnt}.quadrature_ = {1, qx(cnt,:), wx(cnt)};
                    end
                end
                
               this.integral_unit_ = new_integral_unit;  
               this.num_integral_unit_ = total_num_quad_point;
            end
            
            status = true;
        end    
        
        function modifiedQuadPoint(this, integral_unit, nurbs_data)
            import Utility.NurbsUtility.NurbsType
            switch nurbs_data.type_
                case NurbsType.Curve
                    xi = nurbs_data.parametric_mapping_{1};
                    eta = nurbs_data.parametric_mapping_{2};
                    u = integral_unit.quadrature_{2};
                    
                    if length(xi) == 2
                        le = xi(2)-xi(1);
                        xi = xi(1)+ le*u;              
                        integral_unit.quadrature_{2} = [xi eta*ones(size(xi))];
                        integral_unit.quadrature_{3} = integral_unit.quadrature_{3}* le;
                    else
                        le = eta(2)-eta(1);
                        eta = eta(1)+ le*u;              
                        integral_unit.quadrature_{2} = [xi*ones(size(eta)) eta];
                        integral_unit.quadrature_{3} = integral_unit.quadrature_{3}* le;
                    end 
                case NurbsType.Surface
                
                otherwise   
                    disp('input boundary nurbs should be either Curve or Surface!');
            end
        end
        
        
        function generateIntegralUnitByKnotMesh(this,nurbs_data)
            import IntegrationRule.IGA.IntegralUnit
            
            uniqued_knots = nurbs_data.knot_vectors_;
            for i = 1:length(uniqued_knots)
                uniqued_knots{i} = unique(uniqued_knots{i});
            end
            
            dimension = length(uniqued_knots);
            num_element = ones(1, dimension);
            for i = 1:dimension
                num_element(i) = size(uniqued_knots{i},2)-1;
            end
            
            % create integration unit container
            this.num_integral_unit_ = prod(num_element);
            this.integral_unit_ = cell(this.num_integral_unit_, 1);
            
            switch dimension
                case 1
                    for i = 1:num_element(1)
                        element_span = {uniqued_knots{1}(i:i+1)};
                        % new integral unit
                        int_unit = IntegralUnit(this.integral_patch_, element_span);
                        this.integral_unit_{i} = int_unit;
                    end
                case 2
                    for i = 1:num_element(1)
                        for j = 1:num_element(2)
                            n = (i-1)*num_element(2) + j;
                            element_span = {uniqued_knots{1}(i:i+1), uniqued_knots{2}(j:j+1)};
                            % new integral unit
                            int_unit = IntegralUnit(this.integral_patch_, element_span);
                            this.integral_unit_{n} = int_unit;
                        end
                    end
                case 3
                    for i = 1:num_element(1)
                        for j = 1:num_element(2)
                            for k = 1:num_element(3)
                                n = (i-1)*num_element(2)*num_element(3) + (j-1)*num_element(3) + k;
                                element_span = {uniqued_knots{1}(i:i+1), uniqued_knots{2}(j:j+1), uniqued_knots{3}(k:k+1)};
                                % new integral unit
                                int_unit = IntegralUnit(this.integral_patch_, element_span);
                                this.integral_unit_{n} = int_unit;
                            end
                        end
                    end
            end % end switch
        end
        
    end
end

