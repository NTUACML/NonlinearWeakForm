classdef SaintVenantKirchhoff < MaterialBank.MaterialBase 
    properties (Access = private)
        nu_
        E_
        type_
              
        C_
        D_temp_
    end
    
    methods
        function this = SaintVenantKirchhoff(property)
           this@MaterialBank.MaterialBase();
           this.nu_ = property{1};
           this.E_ = property{2};
           this.type_ =  property{3};
           
           this.C_ = [];
           this.D_temp_ = [];
        end
        
        function evaluate(this, F)         
            this.C_ = F'*F;            
        end
        
        function piola_stress = PiolaStress(this)  
            import Utility.BasicUtility.ElasticMaterialType
            
            green_strain = 0.5*(this.C_ - eye(size(this.C_)));
            
            
            if isempty(this.D_temp_)
                this.D_temp_ = this.materialMatrix(); 
            end
           
            switch this.type_
                case ElasticMaterialType.Solid
                    green_strain = [green_strain(1,1), green_strain(2,2), green_strain(3,3), 2*green_strain(1,2), 2*green_strain(2,3), 2*green_strain(1,3)]';
                    piola_stress = this.D_temp_ * green_strain;
                case ElasticMaterialType.PlaneStress
                    green_strain = [green_strain(1,1), green_strain(2,2), 2*green_strain(1,2)]';
                    piola_stress = this.D_temp_ * green_strain;
                case ElasticMaterialType.PlaneStrain
                    green_strain = [green_strain(1,1), green_strain(2,2), 2*green_strain(1,2)]';
                    piola_stress = this.D_temp_ * green_strain;
            end
        end
        
        function D_matrix = materialMatrix(this)
            import Utility.BasicUtility.ElasticMaterialType         
            if isempty(this.D_temp_)
                switch this.type_
                    case ElasticMaterialType.Solid                                  
                        G = this.E_/2/(1+this.nu_);
                        lambda = G/(1-2*this.nu_);
                        
                        D_matrix = zeros(6,6);  
                        
                        D_matrix(1:3, 1:3) = lambda;
                        
                        D_matrix(1,1) = D_matrix(1,1) + 2*G;
                        D_matrix(2,2) = D_matrix(2,2) + 2*G;
                        D_matrix(3,3) = D_matrix(3,3) + 2*G;
                        D_matrix(4,4) = G;
                        D_matrix(5,5) = G;
                        D_matrix(6,6) = G;                       
                    case ElasticMaterialType.PlaneStress
                        temp = this.E_/(1-this.nu_*this.nu_);
                        D_matrix = temp*[1  this.nu_ 0;
                            this.nu_ 1  0;
                            0  0  0.5*(1-this.nu_)];
                    case ElasticMaterialType.PlaneStrain
                        temp = this.E_/(1+this.nu_)/(1-2*this.nu_);
                        D_matrix = temp*[1-this.nu_ this.nu_   0;
                            this.nu_   1-this.nu_ 0;
                            0    0    0.5*(1-2*this.nu_)];
                end
                this.D_temp_ = D_matrix;
            else
                D_matrix = this.D_temp_;
            end
        end
    end
    
    
    
    
    
    
end

