function XML_Writer
clear; clc; home

import Utility.Resources.XML_DataMaker.IGA.NurbsGenerator
import Utility.Resources.NurbsToolBox.*

%% Generate nurbs object
% case_name = 'Unit_Square';
% case_name = 'Rectangle';
% case_name = 'Plane_quarter_hole';
% case_name = 'Plane_quarter_hole_2';
% case_name = 'Solid_quarter_hole';
% case_name = 'Lens_top_right'; 
% case_name = 'Lens_bottom_right'; 
% case_name = 'Lens_top_left'; 
% case_name = 'Lens_bottom_left'; 
% case_name = '3D_Lens_left'; 
% case_name = '3D_Lens_right';
case_name = 'Unit_Cube';

nurbs = NurbsGenerator(case_name);

%% Document setting
project_name = 'ArtPDE';
doc_format = ['IGA_', case_name];
version = '1.1';

import Utility.NurbsUtility.NurbsType
switch nurbs.type_
    case NurbsType.Solid
        geo_dim = '3';
    case NurbsType.Plane
        geo_dim = '2';
    case NurbsType.Surface
        geo_dim = '2';
    case NurbsType.Curve
        geo_dim = '1';
end

%% Document create
doc_geo = com.mathworks.xml.XMLUtils.createDocument('Geometry');
% doc_init = com.mathworks.xml.XMLUtils.createDocument('Initial');
% doc_mat = com.mathworks.xml.XMLUtils.createDocument('Material');

%% Document node
doc_geo_node = doc_geo.getDocumentElement();
doc_geo_node.setAttribute('version',version);
doc_geo_node.setAttribute('dim', geo_dim);

% doc_init_node = doc_init.getDocumentElement();
% doc_init_node.setAttribute('version',version);
% 
% doc_mat_node = doc_mat.getDocumentElement();
% doc_mat_node.setAttribute('version',version);

%% Data create (Geometry part)
doc_handle = doc_geo;

% /Unit
unit_node = DataNodeCreate('Unit', doc_geo_node, doc_handle);
unit_node.setAttribute('format','NURBS');

% /Unit/Patch (Domain)
patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
patch_node.setAttribute('region','Domain');
patch_node.setAttribute('name','domain_1');
% patch_node.setAttribute('type','Surface');

% /Unit/Patch/ControlPoint
CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
CTPT_node.setAttribute('dim','3');

% /Unit/Patch/ControlPoint/Point
point_data = nurbs.control_points_(:,:);

for i = 1 : size(point_data, 1)
    point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
    point_row_str = num2str(point_data(i, :));
    point_node.appendChild(doc_handle.createTextNode(point_row_str));
end

% /Unit/Patch/Order
order_node = DataNodeCreate('Order', patch_node, doc_handle);
order_str = num2str(nurbs.order_);
order_node.appendChild(doc_handle.createTextNode(order_str));

% /Unit/Patch/Knot
knot_node = DataNodeCreate('Knot', patch_node, doc_handle);

% /Unit/Patch/Knot/Vector
knot_data = nurbs.knot_vectors_;

for i = 1 : size(knot_data, 2)
    vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
    knot_row_str = num2str(knot_data{i});
    vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
    vector_node.setAttribute('dof',num2str(length(knot_data{i})));
end

% /Unit/Patch (Boundary)
switch geo_dim
    case '3'
        import Utility.BasicUtility.TensorProduct       
        
        % Get internal knots
        bool = abs(nurbs.knot_vectors_{1}) < eps | abs(nurbs.knot_vectors_{1}-1) < eps;
        xi_i = nurbs.knot_vectors_{1}(~bool);
        
        bool = abs(nurbs.knot_vectors_{2}) < eps | abs(nurbs.knot_vectors_{2}-1) < eps;
        eta_i = nurbs.knot_vectors_{2}(~bool);
        
        bool = abs(nurbs.knot_vectors_{3}) < eps | abs(nurbs.knot_vectors_{3}-1) < eps;
        zeta_i = nurbs.knot_vectors_{3}(~bool);
        
        % Get order to be elevated
        order_elevated = nurbs.order_-1;
        
        %% CREATE BOUNDARY NURBS PATCH (ZETA = 1)
        patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
        patch_node.setAttribute('region','Boundary');
        patch_node.setAttribute('name','zeta_1');
        
        % /Unit/Patch/ControlPoint
        CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
        CTPT_node.setAttribute('dim','3');

        % /Unit/Patch/ControlPoint/Point
        surf = nrb4surf([0.0 0.0 1.0],[1.0 0.0 1.0],[0.0 1.0 1.0],[1.0 1.0 1.0]);      
        surf = nrbdegelev(surf, [order_elevated(1) order_elevated(2)]);
        surf = nrbkntins(surf,{xi_i eta_i});
      
        TD = TensorProduct({surf.number(1) surf.number(2)});
        for id = 1 : TD.total_num_
            local = TD.to_local_index(id);
            point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
            point_row_str = num2str(surf.coefs(:,local{1}, local{2})');
            point_node.appendChild(doc_handle.createTextNode(point_row_str));
        end
        
        % /Unit/Patch/Order
        order_node = DataNodeCreate('Order', patch_node, doc_handle);
        order_str = num2str(surf.order-1);
        order_node.appendChild(doc_handle.createTextNode(order_str));
        
        % /Unit/Patch/Knot
        knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
        
        % /Unit/Patch/Knot/Vector
        knot_data = surf.knots;
        for i = 1 : size(knot_data, 2)
            vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
            knot_row_str = num2str(knot_data{i});
            vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
            vector_node.setAttribute('dof',num2str(length(knot_data{i})));
        end 
        
        %% CREATE BOUNDARY NURBS PATCH (ZETA = 0)
        patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
        patch_node.setAttribute('region','Boundary');
        patch_node.setAttribute('name','zeta_0');
        
        % /Unit/Patch/ControlPoint
        CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
        CTPT_node.setAttribute('dim','3');

        % /Unit/Patch/ControlPoint/Point
        surf = nrb4surf([0.0 0.0 0.0],[1.0 0.0 0.0],[0.0 1.0 0.0],[1.0 1.0 0.0]);
        surf = nrbdegelev(surf, [order_elevated(1) order_elevated(2)]);
        surf = nrbkntins(surf,{xi_i eta_i});
      
        TD = TensorProduct({surf.number(1) surf.number(2)});
        for id = 1 : TD.total_num_
            local = TD.to_local_index(id);
            point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
            point_row_str = num2str(surf.coefs(:,local{1}, local{2})');
            point_node.appendChild(doc_handle.createTextNode(point_row_str));
        end
        
        % /Unit/Patch/Order
        order_node = DataNodeCreate('Order', patch_node, doc_handle);
        order_str = num2str(surf.order-1);
        order_node.appendChild(doc_handle.createTextNode(order_str));
        
        % /Unit/Patch/Knot
        knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
        
        % /Unit/Patch/Knot/Vector
        knot_data = surf.knots;
        for i = 1 : size(knot_data, 2)
            vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
            knot_row_str = num2str(knot_data{i});
            vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
            vector_node.setAttribute('dof',num2str(length(knot_data{i})));
        end 
        
        %% CREATE BOUNDARY NURBS PATCH (ETA = 1)
        patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
        patch_node.setAttribute('region','Boundary');
        patch_node.setAttribute('name','eta_1');
        
        % /Unit/Patch/ControlPoint
        CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
        CTPT_node.setAttribute('dim','3');

        % /Unit/Patch/ControlPoint/Point
        surf = nrb4surf([0.0 1.0 0.0],[0.0 1.0 1.0],[1.0 1.0 0.0],[1.0 1.0 1.0]);
        surf = nrbdegelev(surf, [order_elevated(1) order_elevated(3)]);
        surf = nrbkntins(surf,{xi_i zeta_i});
        
        TD = TensorProduct({surf.number(1) surf.number(2)});
        for id = 1 : TD.total_num_
            local = TD.to_local_index(id);
            point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
            point_row_str = num2str(surf.coefs(:,local{1}, local{2})');
            point_node.appendChild(doc_handle.createTextNode(point_row_str));
        end
        
        % /Unit/Patch/Order
        order_node = DataNodeCreate('Order', patch_node, doc_handle);
        order_str = num2str(surf.order-1);
        order_node.appendChild(doc_handle.createTextNode(order_str));
        
        % /Unit/Patch/Knot
        knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
        
        % /Unit/Patch/Knot/Vector
        knot_data = surf.knots;
        for i = 1 : size(knot_data, 2)
            vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
            knot_row_str = num2str(knot_data{i});
            vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
            vector_node.setAttribute('dof',num2str(length(knot_data{i})));
        end 
        
        %% CREATE BOUNDARY NURBS PATCH (ETA = 0)
        patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
        patch_node.setAttribute('region','Boundary');
        patch_node.setAttribute('name','eta_0');
        
        % /Unit/Patch/ControlPoint
        CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
        CTPT_node.setAttribute('dim','3');

        % /Unit/Patch/ControlPoint/Point
        surf = nrb4surf([0.0 0.0 0.0],[0.0 0.0 1.0],[1.0 0.0 0.0],[1.0 0.0 1.0]);
        surf = nrbdegelev(surf, [order_elevated(1) order_elevated(3)]);
        surf = nrbkntins(surf,{xi_i zeta_i});
      
        TD = TensorProduct({surf.number(1) surf.number(2)});
        for id = 1 : TD.total_num_
            local = TD.to_local_index(id);
            point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
            point_row_str = num2str(surf.coefs(:,local{1}, local{2})');
            point_node.appendChild(doc_handle.createTextNode(point_row_str));
        end
        
        % /Unit/Patch/Order
        order_node = DataNodeCreate('Order', patch_node, doc_handle);
        order_str = num2str(surf.order-1);
        order_node.appendChild(doc_handle.createTextNode(order_str));
        
        % /Unit/Patch/Knot
        knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
        
        % /Unit/Patch/Knot/Vector
        knot_data = surf.knots;
        for i = 1 : size(knot_data, 2)
            vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
            knot_row_str = num2str(knot_data{i});
            vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
            vector_node.setAttribute('dof',num2str(length(knot_data{i})));
        end 
        
        %% CREATE BOUNDARY NURBS PATCH (XI = 1)
        patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
        patch_node.setAttribute('region','Boundary');
        patch_node.setAttribute('name','xi_1');
        
        % /Unit/Patch/ControlPoint
        CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
        CTPT_node.setAttribute('dim','3');

        % /Unit/Patch/ControlPoint/Point
        surf = nrb4surf([1.0 0.0 0.0],[1.0 1.0 0.0],[1.0 0.0 1.0],[1.0 1.0 1.0]);
        surf = nrbdegelev(surf, [order_elevated(2) order_elevated(3)]);
        surf = nrbkntins(surf,{eta_i zeta_i});
      
        TD = TensorProduct({surf.number(1) surf.number(2)});
        for id = 1 : TD.total_num_
            local = TD.to_local_index(id);
            point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
            point_row_str = num2str(surf.coefs(:,local{1}, local{2})');
            point_node.appendChild(doc_handle.createTextNode(point_row_str));
        end
        
        % /Unit/Patch/Order
        order_node = DataNodeCreate('Order', patch_node, doc_handle);
        order_str = num2str(surf.order-1);
        order_node.appendChild(doc_handle.createTextNode(order_str));
        
        % /Unit/Patch/Knot
        knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
        
        % /Unit/Patch/Knot/Vector
        knot_data = surf.knots;
        for i = 1 : size(knot_data, 2)
            vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
            knot_row_str = num2str(knot_data{i});
            vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
            vector_node.setAttribute('dof',num2str(length(knot_data{i})));
        end 
        
        %% CREATE BOUNDARY NURBS PATCH (XI = 0)
        patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
        patch_node.setAttribute('region','Boundary');
        patch_node.setAttribute('name','xi_0');
        
        % /Unit/Patch/ControlPoint
        CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
        CTPT_node.setAttribute('dim','3');

        % /Unit/Patch/ControlPoint/Point
        surf = nrb4surf([0.0 0.0 0.0],[0.0 1.0 0.0],[0.0 0.0 1.0],[0.0 1.0 1.0]);
        surf = nrbdegelev(surf, [order_elevated(2) order_elevated(3)]);
        surf = nrbkntins(surf,{eta_i zeta_i});
      
        TD = TensorProduct({surf.number(1) surf.number(2)});
        for id = 1 : TD.total_num_
            local = TD.to_local_index(id);
            point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
            point_row_str = num2str(surf.coefs(:,local{1}, local{2})');
            point_node.appendChild(doc_handle.createTextNode(point_row_str));
        end
        
        % /Unit/Patch/Order
        order_node = DataNodeCreate('Order', patch_node, doc_handle);
        order_str = num2str(surf.order-1);
        order_node.appendChild(doc_handle.createTextNode(order_str));
        
        % /Unit/Patch/Knot
        knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
        
        % /Unit/Patch/Knot/Vector
        knot_data = surf.knots;
        for i = 1 : size(knot_data, 2)
            vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
            knot_row_str = num2str(knot_data{i});
            vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
            vector_node.setAttribute('dof',num2str(length(knot_data{i})));
        end 
        
    case '2'      
        % Get internal knots 
        bool = abs(nurbs.knot_vectors_{1}) < eps | abs(nurbs.knot_vectors_{1}-1) < eps;
        xi_i = nurbs.knot_vectors_{1}(~bool);
        
        bool = abs(nurbs.knot_vectors_{2}) < eps | abs(nurbs.knot_vectors_{2}-1) < eps;
        eta_i = nurbs.knot_vectors_{2}(~bool);
 
        % Get order to be elevated
        order_elevated = nurbs.order_-1;
               
        % CREATE BOUNDARY NURBS PATCH (ETA = 1)
        patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
        patch_node.setAttribute('region','Boundary');
        patch_node.setAttribute('name','eta_1');
        
        % /Unit/Patch/ControlPoint
        CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
        CTPT_node.setAttribute('dim','3');

        % /Unit/Patch/ControlPoint/Point
        crv = nrbline([0 1 0]', [1 1 0]');
        crv = nrbdegelev(crv, order_elevated(1));
        crv = nrbkntins(crv, xi_i);
        
        point_data = crv.coefs';       
        for i = 1 : size(point_data, 1)
            point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
            point_row_str = num2str(point_data(i, :));
            point_node.appendChild(doc_handle.createTextNode(point_row_str));
        end
        
        % /Unit/Patch/Order
        order_node = DataNodeCreate('Order', patch_node, doc_handle);
        order_str = num2str(crv.order-1);
        order_node.appendChild(doc_handle.createTextNode(order_str));
        
        % /Unit/Patch/Knot
        knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
        
        % /Unit/Patch/Knot/Vector
        knot_data = {crv.knots};
        for i = 1 : size(knot_data, 1)
            vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
            knot_row_str = num2str(knot_data{i});
            vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
            vector_node.setAttribute('dof',num2str(length(knot_data{i})));
        end 
        
        % CREATE BOUNDARY NURBS PATCH (ETA = 0)
        patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
        patch_node.setAttribute('region','Boundary');
        patch_node.setAttribute('name','eta_0');
        
        % /Unit/Patch/ControlPoint
        CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
        CTPT_node.setAttribute('dim','3');

        % /Unit/Patch/ControlPoint/Point
        crv = nrbline([0 0 0]',[1 0 0]');
        crv = nrbdegelev(crv, order_elevated(1));
        crv = nrbkntins(crv, xi_i);
      
        point_data = crv.coefs';                            
        for i = 1 : size(point_data, 1)
            point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
            point_row_str = num2str(point_data(i, :));
            point_node.appendChild(doc_handle.createTextNode(point_row_str));
        end
        
        % /Unit/Patch/Order
        order_node = DataNodeCreate('Order', patch_node, doc_handle);
        order_str = num2str(crv.order-1);
        order_node.appendChild(doc_handle.createTextNode(order_str));
        
        % /Unit/Patch/Knot
        knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
        
        % /Unit/Patch/Knot/Vector
        knot_data = {crv.knots};
        for i = 1 : size(knot_data, 1)
            vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
            knot_row_str = num2str(knot_data{i});
            vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
            vector_node.setAttribute('dof',num2str(length(knot_data{i})));
        end 
        
        % CREATE BOUNDARY NURBS PATCH (XI = 0)
        patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
        patch_node.setAttribute('region','Boundary');
        patch_node.setAttribute('name','xi_0');
        
        % /Unit/Patch/ControlPoint
        CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
        CTPT_node.setAttribute('dim','3');

        % /Unit/Patch/ControlPoint/Point
        crv = nrbline([0 0 0]',[0 1 0]');
        crv = nrbdegelev(crv, order_elevated(2));
        crv = nrbkntins(crv, eta_i);
      
        point_data = crv.coefs';                   
        for i = 1 : size(point_data, 1)
            point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
            point_row_str = num2str(point_data(i, :));
            point_node.appendChild(doc_handle.createTextNode(point_row_str));
        end
        
        % /Unit/Patch/Order
        order_node = DataNodeCreate('Order', patch_node, doc_handle);
        order_str = num2str(crv.order-1);
        order_node.appendChild(doc_handle.createTextNode(order_str));
        
        % /Unit/Patch/Knot
        knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
        
        % /Unit/Patch/Knot/Vector
        knot_data = {crv.knots};
        for i = 1 : size(knot_data, 1)
            vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
            knot_row_str = num2str(knot_data{i});
            vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
            vector_node.setAttribute('dof',num2str(length(knot_data{i})));
        end 
        
        % CREATE RIGHT BOUNDARY NURBS PATCH (XI = 1)
        patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
        patch_node.setAttribute('region','Boundary');
        patch_node.setAttribute('name','xi_1');
        
        % /Unit/Patch/ControlPoint
        CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
        CTPT_node.setAttribute('dim','3');

        % /Unit/Patch/ControlPoint/Point
        crv = nrbline([1 0 0]',[1 1 0]');
        crv = nrbdegelev(crv, order_elevated(2));
        crv = nrbkntins(crv, eta_i);
      
        point_data = crv.coefs';
        for i = 1 : size(point_data, 1)
            point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
            point_row_str = num2str(point_data(i, :));
            point_node.appendChild(doc_handle.createTextNode(point_row_str));
        end
        
        % /Unit/Patch/Order
        order_node = DataNodeCreate('Order', patch_node, doc_handle);
        order_str = num2str(crv.order-1);
        order_node.appendChild(doc_handle.createTextNode(order_str));
        
        % /Unit/Patch/Knot
        knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
        
        % /Unit/Patch/Knot/Vector
        knot_data = {crv.knots};
        for i = 1 : size(knot_data, 1)
            vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
            knot_row_str = num2str(knot_data{i});
            vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
            vector_node.setAttribute('dof',num2str(length(knot_data{i})));
        end 
        
        
%         % START Test Trimmed Curve
%         % CREATE TRIMMED BOUNDARY NURBS PATCH
%         patch_node = DataNodeCreate('Patch', unit_node, doc_handle);
%         patch_node.setAttribute('region','Boundary');
%         patch_node.setAttribute('name','trimmed_curve');
%         
%         % /Unit/Patch/ControlPoint
%         CTPT_node = DataNodeCreate('ControlPoint', patch_node, doc_handle);
%         CTPT_node.setAttribute('dim','3');
% 
%         % /Unit/Patch/ControlPoint/Point
%         pnts = [0.1 0.3 0.6 1.0; 
%                 0.0 0.5 0.6 0.8; 
%                 0.0 0.0 0.0 0.0]; 
%             
%         crv = nrbmak(pnts,[0 0 0 1/2 1 1 1]); 
%         
%         point_data = crv.coefs';
%         for i = 1 : size(point_data, 1)
%             point_node = DataNodeCreate('Point', CTPT_node, doc_handle);
%             point_row_str = num2str(point_data(i, :));
%             point_node.appendChild(doc_handle.createTextNode(point_row_str));
%         end
%         
%         % /Unit/Patch/Order
%         order_node = DataNodeCreate('Order', patch_node, doc_handle);
%         order_str = num2str(crv.order-1);
%         order_node.appendChild(doc_handle.createTextNode(order_str));
%         
%         % /Unit/Patch/Knot
%         knot_node = DataNodeCreate('Knot', patch_node, doc_handle);
%         
%         % /Unit/Patch/Knot/Vector
%         knot_data = {crv.knots};
%         for i = 1 : size(knot_data, 1)
%             vector_node = DataNodeCreate('Vector', knot_node, doc_handle);
%             knot_row_str = num2str(knot_data{i});
%             vector_node.appendChild(doc_handle.createTextNode(knot_row_str));
%             vector_node.setAttribute('dof',num2str(length(knot_data{i})));
%         end
%         % END Test Trimmed Curve
        
        
        
end

%% Document write
file_name_geo = [project_name, '_', doc_format ,'.art_geometry'];
xmlwrite(file_name_geo, doc_geo);

% file_name_init = [project_name, '_', doc_format,'.art_initial'];
% xmlwrite(file_name_init, doc_init);
% 
% file_name_mat = [project_name, '_', doc_format,'.art_material'];
% xmlwrite(file_name_mat, doc_mat);
end

%% Using Function
function current_node = DataNodeCreate(name, upper_node, doc_handle)
    current_node = doc_handle.createElement(name);
    upper_node.appendChild(current_node);
end