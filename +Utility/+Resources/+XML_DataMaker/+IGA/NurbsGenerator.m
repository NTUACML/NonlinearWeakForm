% Function used to generate the nurbs geometry.
% You ca modify the knot or order by using the following utility functions:
% 1. nurbs = degreeElevation(nurbs, [1 1]);
% 2. nurbs = knotInsertion(nurbs, {[0.2 0.8] [0.3 0.7]});

function nurbs_out = NurbsGenerator(nurbs_name)
import Utility.NurbsUtility.Nurbs

switch nurbs_name
    case 'Unit_Square'
        nurbs = Rectangle(1, 1, [0.5 0.5]);
%         nurbs = degreeElevation(nurbs, [2 2]);
%         n = 10;
%         t = linspace(1/n, 1-1/n, n-1);
%         nurbs = knotInsertion(nurbs, {t t});
        status = true;
    case 'Rectangle'
        % Uniaxial tensor beam
%         D = 0.5; L = 1;
%         nurbs = Rectangle(D, L, [L/2 D/2]);       
%         nurbs = degreeElevation(nurbs, [1 1]);
%         n = 5;
%         t_1 = linspace(1/n, 1-1/n, n-1);
%         n = 2;
%         t_2 = linspace(1/n, 1-1/n, n-1);

        % Canteliver beam
        D = 0.2; L = 1;
        nurbs = Rectangle(D, L, [L/2 0]);
        nurbs = degreeElevation(nurbs, [1 1]);
        n = 5;
        t_1 = linspace(1/n, 1-1/n, n-1);
%         n = 2;
%         t_2 = linspace(1/n, 1-1/n, n-1); 
        nurbs = knotInsertion(nurbs, {t_1 [0.5 0.5]});
        status = true;  
    case 'Plane_quarter_hole'
        nurbs = Plane_quarter_hole();
        nurbs = degreeElevation(nurbs, [1 1]);
        
        nurbs = knotInsertion(nurbs, {[0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9], 0.1:0.1:0.9});     
        status = true;
    case 'Plane_quarter_hole_2'
        nurbs = Plane_quarter_hole_2();
%         nurbs = degreeElevation(nurbs, [1 1]);
        nurbs = knotInsertion(nurbs, {[0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9], 0.1:0.1:0.9}); 
        status = true;    
    case 'Solid_quarter_hole'        
        nurbs = Solid_quater_hole(1);
        status = true;    
    case 'Lens_top_right'
        nurbs = Surface_Lens('top_right');
        status = true;
    case 'Lens_bottom_right'
        nurbs = Surface_Lens('bottom_right');
        status = true;
    case 'Lens_top_left'
        nurbs = Surface_Lens('top_left');
        status = true;
    case 'Lens_bottom_left'
        nurbs = Surface_Lens('bottom_left');
        status = true;    
    case '3D_Lens_left'
        nurbs = Solid_Lens('left');
        status = true;
    case '3D_Lens_right'
        nurbs = Solid_Lens('right');
        status = true; 
    case 'Unit_Cube'
%         nurbs = Solid_Cube(1, 1, 1, [0.5 0.5 0.5]);
%         nurbs = degreeElevation(nurbs, [1 1 1]);
%         n = 2;
%         t = linspace(1/n, 1-1/n, n-1);      
%         nurbs = knotInsertion(nurbs, {t t t});
        nurbs = Solid_Cube(6, 0.1, 0.2, 0.5*[6, 0.1, 0.2]);
%         nurbs = degreeElevation(nurbs, [1 1 1]);
%         n = 10;
%         t = linspace(1/n, 1-1/n, n-1);
%         nurbs = knotInsertion(nurbs, {t [] []});
        status = true;        
    otherwise
        str = [nurbs_name, 'does not exist in the current library.'];
        disp(str);
        status = false;
end

if status
    if isa(nurbs, 'Nurbs')
        nurbs_out = nurbs;
    else
        nurbs_out = OutPutFiles(nurbs);
    end
end

end

function nurbs_out= OutPutFiles(nurbs)
import Utility.BasicUtility.TensorProduct

% nurbs knots
knots = nurbs.knots;

% nurbs control point
geo_dim = size(knots, 2);
control_pnt = zeros(prod(nurbs.number), 4);
switch geo_dim
    case 1
        TD_1 = TensorProduct({nurbs.number(1)});
        for k = 1:TD_1.total_num_
            local_index = TD_1.to_local_index(k);
            control_pnt(k,:) = nurbs.coefs(:, local_index{1})';
        end
    case 2
        TD_2 = TensorProduct({nurbs.number(1) nurbs.number(2)});
        for k = 1:TD_2.total_num_
            local_index = TD_2.to_local_index(k);
            control_pnt(k,:) = nurbs.coefs(:, local_index{1}, local_index{2})';
        end
    case 3
        TD_3 = TensorProduct({nurbs.number(1) nurbs.number(2) nurbs.number(3)});
        for k = 1:TD_3.total_num_
            local_index = TD_3.to_local_index(k);
            control_pnt(k,:) = nurbs.coefs(:, local_index{1}, local_index{2}, local_index{3})';
        end
end

control_pnt(:,1) = control_pnt(:,1)./control_pnt(:,4);
control_pnt(:,2) = control_pnt(:,2)./control_pnt(:,4);
control_pnt(:,3) = control_pnt(:,3)./control_pnt(:,4);

% nurbs order
order = nurbs.order-1;

import Utility.BasicUtility.PointList
control_point_list = PointList(control_pnt);

import Utility.NurbsUtility.Nurbs
nurbs_out = Nurbs(knots, order, control_point_list);

end

% Supported geometries
function nurbs = Plane_quarter_hole
import Utility.BasicUtility.PointList
import Utility.NurbsUtility.Nurbs

knot_vectors = {[0 0 0 0.5 1 1 1]  [0 0 0 1 1 1]};
order = [2 2];

temp = [ 
    0 1 0 1;
    -1+sqrt(2) 1 0 (1+1/sqrt(2))/2;
    1 sqrt(2)-1 0 (1+1/sqrt(2))/2;
    1 0 0 1;
 
    0 2.5 0 1;
    0.75 2.5 0 1;
    2.5 0.75 0 1;
    2.5 0 0 1;
 
    0 4 0 1;
    4 4 0 1;
    4 4 0 1;
    4 0 0 1
    ];

control_point_list = PointList(temp);

nurbs = Nurbs(knot_vectors, order, control_point_list);

end

function nurbs = Plane_quarter_hole_2
import Utility.BasicUtility.PointList
import Utility.NurbsUtility.Nurbs

knot_vectors = {[0 0 0 0.5 1 1 1]  [0 0 0 1 1 1]};
order = [2 2];

temp = [ 
    0 1 0 1;
    -1+sqrt(2) 1 0 (1+1/sqrt(2))/2;
    1 sqrt(2)-1 0 (1+1/sqrt(2))/2;
    1 0 0 1;
 
    0 2.5 0 1;
    0.75 2.5 0 1;
    2.5 0.75 0 1;
    2.5 0 0 1;
 
    0 4 0 1;
    4 4 0 1;
    4 4 0 1;
    4 0 0 1
    ];


for i = 1:size(temp,1)
    temp(i,1:2) = temp(i,1:2)*[0 1; -1 0];
end

control_point_list = PointList(temp);

nurbs = Nurbs(knot_vectors, order, control_point_list);

end

function nurbs = Solid_quater_hole(height)
import Utility.BasicUtility.PointList
import Utility.NurbsUtility.Nurbs
import Utility.NurbsUtility.NurbsType
nurbs = Plane_quarter_hole();

nurbs.basis_number_(3) = 2;
nurbs.knot_vectors_{1,3} = [0 0 1 1];
nurbs.order_(3) = 1;

temp = nurbs.control_points_(:,:);
temp(:,3) = height;
nurbs.control_points_ = PointList([nurbs.control_points_(:,:); temp]);

nurbs.type_ = NurbsType.Solid;
end

function nurbs = Rectangle(width, length, center)
import Utility.BasicUtility.PointList
import Utility.NurbsUtility.Nurbs

knot_vectors = {[0 0 1 1] [0 0 1 1]};
order = [1 1];
temp = [center(1)-length/2 center(2)-width/2 0.0 1.0;
        center(1)+length/2 center(2)-width/2 0.0 1.0;
        center(1)-length/2 center(2)+width/2 0.0 1.0;
        center(1)+length/2 center(2)+width/2 0.0 1.0];
control_point_list = PointList(temp);

nurbs = Nurbs(knot_vectors, order, control_point_list);
end

function nurbs = Solid_Lens(str)
import Utility.BasicUtility.PointList
import Utility.NurbsUtility.Nurbs

knot_vectors = {[0 0 0 1 1 1] [0 0 0 0.5 0.5 1 1 1] [0 0 0 0 1 1 1 1]};
order = [2 2 3];

switch str
    case 'right'
        temp = [89.673928030621 -0.000000000000 -6.879260917442 1.000000000000;
            95.904289874264 -0.000000000000 -8.010000000000 0.995973631625;
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000;
            89.673928030621 12.562500469379 -6.879260917442 0.707106781187;
            95.904289874264 6.332138625736 -8.010000000000 0.704259708805;
            102.236428500000 -0.000000000000 -8.009999999999 0.707106781187;
            102.236428500000 12.562500469379 -6.879260917442 1.000000000000;
            102.236428500000 6.332138625736 -8.010000000000 0.995973631625;
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000;
            114.798928969379 12.562500469379 -6.879260917442 0.707106781187;
            108.568567125736 6.332138625736 -8.010000000000 0.704259708805;
            102.236428500000 -0.000000000000 -8.009999999999 0.707106781187;
            114.798928969379 -0.000000000000 -6.879260917442 1.000000000000;
            108.568567125736 -0.000000000000 -8.010000000000 0.995973631625;
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000;
            
            89.673928030621	-0.000000000000	-6.586173945	1.000000000000;
            95.904289874264	-0.000000000000	-7.34	0.995973631625;
            102.236428500000	-0.000000000000	-7.34	1.000000000000;
            89.673928030621	12.562500469379	-6.586173945	0.707106781187;
            95.904289874264	6.332138625736	-7.34	0.704259708805;
            102.236428500000	-0.000000000000	-7.34	0.707106781187;
            102.236428500000	12.562500469379	-6.586173945	1.000000000000;
            102.236428500000	6.332138625736	-7.34	0.995973631625;
            102.236428500000	-0.000000000000	-7.34	1.000000000000;
            114.798928969379	12.562500469379	-6.586173945	0.707106781187;
            108.568567125736	6.332138625736	-7.34	0.704259708805;
            102.236428500000	-0.000000000000	-7.34	0.707106781187;
            114.798928969379	-0.000000000000	-6.586173945	1.000000000000;
            108.568567125736	-0.000000000000	-7.34	0.995973631625;
            102.236428500000	-0.000000000000	-7.34	1.000000000000;
            
            89.673928030621	-0.000000000000	-6.293086972	1.000000000000;
            95.904289874264	-0.000000000000	-6.67	0.995973631625;
            102.236428500000	-0.000000000000	-6.67	1.000000000000;
            89.673928030621	12.562500469379	-6.293086972	0.707106781187;
            95.904289874264	6.332138625736	-6.67	0.704259708805;
            102.236428500000	-0.000000000000	-6.67	0.707106781187;
            102.236428500000	12.562500469379	-6.293086972	1.000000000000;
            102.236428500000	6.332138625736	-6.67	0.995973631625
            102.236428500000	-0.000000000000	-6.67	1.000000000000
            114.798928969379	12.562500469379	-6.293086972	0.707106781187
            108.568567125736	6.332138625736	-6.67	0.704259708805;
            102.236428500000	-0.000000000000	-6.67	0.707106781187;
            114.798928969379	-0.000000000000	-6.293086972	1.000000000000;
            108.568567125736	-0.000000000000	-6.67	0.995973631625;
            102.236428500000	-0.000000000000	-6.67	1.000000000000;
            
            89.673928030621	-0.000000000000	-6.000000000000	1.000000000000;
            95.904289874264	-0.000000000000	-6.000000000000	0.995973631625;
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000;
            89.673928030621	12.562500469379	-6.000000000000	0.707106781187;
            95.904289874264	6.332138625736	-6.000000000000	0.704259708805;
            102.236428500000	-0.000000000000	-6.000000000000	0.707106781187;
            102.236428500000	12.562500469379	-6.000000000000	1.000000000000;
            102.236428500000	6.332138625736	-6.000000000000	0.995973631625;
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000;
            114.798928969379	12.562500469379	-6.000000000000	0.707106781187;
            108.568567125736	6.332138625736	-6.000000000000	0.704259708805;
            102.236428500000	-0.000000000000	-6.000000000000	0.707106781187;
            114.798928969379	-0.000000000000	-6.000000000000	1.000000000000;
            108.568567125736	-0.000000000000	-6.000000000000	0.995973631625;
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000];
    case 'left'
        temp = [
            114.798928969379 -0.000000000000 -6.879260917442 1.000000000000
            108.568567125736 -0.000000000000 -8.010000000000 0.995973631625
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000
            114.798928969379 -12.562500469379 -6.879260917442 0.707106781187
            108.568567125736 -6.332138625736 -8.010000000000 0.704259708805
            102.236428500000 -0.000000000000 -8.009999999999 0.707106781187
            102.236428500000 -12.562500469379 -6.879260917442 1.000000000000
            102.236428500000 -6.332138625736 -8.010000000000 0.995973631625
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000
            89.673928030621 -12.562500469379 -6.879260917442 0.707106781187
            95.904289874264 -6.332138625736 -8.010000000000 0.704259708805
            102.236428500000 -0.000000000000 -8.009999999999 0.707106781187
            89.673928030621 -0.000000000000 -6.879260917442 1.000000000000
            95.904289874264 -0.000000000000 -8.010000000000 0.995973631625
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000
            
            114.798928969379	-0.000000000000	-6.586173945	1.000000000000
            108.568567125736	-0.000000000000	-7.34	0.995973631625
            102.236428500000	-0.000000000000	-7.34	1.000000000000
            114.798928969379	-12.562500469379	-6.586173945	0.707106781187
            108.568567125736	-6.332138625736	-7.34	0.704259708805
            102.236428500000	-0.000000000000	-7.34	0.707106781187
            102.236428500000	-12.562500469379	-6.586173945	1.000000000000
            102.236428500000	-6.332138625736	-7.34	0.995973631625
            102.236428500000	-0.000000000000	-7.34	1.000000000000
            89.673928030621	-12.562500469379	-6.586173945	0.707106781187
            95.904289874264	-6.332138625736	-7.34	0.704259708805
            102.236428500000	-0.000000000000	-7.34	0.707106781187
            89.673928030621	-0.000000000000	-6.586173945	1.000000000000
            95.904289874264	-0.000000000000	-7.34	0.995973631625
            102.236428500000	-0.000000000000	-7.34	1.000000000000
            
            114.798928969379	-0.000000000000	-6.293086972	1.000000000000
            108.568567125736	-0.000000000000	-6.67	0.995973631625
            102.236428500000	-0.000000000000	-6.67	1.000000000000
            114.798928969379	-12.562500469379	-6.293086972	0.707106781187
            108.568567125736	-6.332138625736	-6.67	0.704259708805
            102.236428500000	-0.000000000000	-6.67	0.707106781187
            102.236428500000	-12.562500469379	-6.293086972	1.000000000000
            102.236428500000	-6.332138625736	-6.67	0.995973631625
            102.236428500000	-0.000000000000	-6.67	1.000000000000
            89.673928030621	-12.562500469379	-6.293086972	0.707106781187
            95.904289874264	-6.332138625736	-6.67	0.704259708805
            102.236428500000	-0.000000000000	-6.67	0.707106781187
            89.673928030621	-0.000000000000	-6.293086972	1.000000000000
            95.904289874264	-0.000000000000	-6.67	0.995973631625
            102.236428500000	-0.000000000000	-6.67	1.000000000000
            
            114.798928969379	-0.000000000000	-6.000000000000	1.000000000000
            108.568567125736	-0.000000000000	-6.000000000000	0.995973631625
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000
            114.798928969379	-12.562500469379	-6.000000000000	0.707106781187
            108.568567125736	-6.332138625736	-6.000000000000	0.704259708805
            102.236428500000	-0.000000000000	-6.000000000000	0.707106781187
            102.236428500000	-12.562500469379	-6.000000000000	1.000000000000
            102.236428500000	-6.332138625736	-6.000000000000	0.995973631625
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000
            89.673928030621	-12.562500469379	-6.000000000000	0.707106781187
            95.904289874264	-6.332138625736	-6.000000000000	0.704259708805
            102.236428500000	-0.000000000000	-6.000000000000	0.707106781187
            89.673928030621	-0.000000000000	-6.000000000000	1.000000000000
            95.904289874264	-0.000000000000	-6.000000000000	0.995973631625
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000];
end
control_point_list = PointList(temp);

nurbs = Nurbs(knot_vectors, order, control_point_list);
end

function nurbs = Surface_Lens(str)
import Utility.BasicUtility.PointList
import Utility.NurbsUtility.Nurbs

knot_vectors = {[0 0 0 1 1 1] [0 0 0 0.5 0.5 1 1 1]};
order = [2 2];

switch str
    case 'bottom_right'
        temp = [114.798928969379	-0.000000000000	-6.000000000000	1.000000000000
            108.568567125736	-0.000000000000	-6.000000000000	0.995973631625
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000
            114.798928969379	-12.562500469379	-6.000000000000	0.707106781187
            108.568567125736	-6.332138625736	-6.000000000000	0.704259708805
            102.236428500000	-0.000000000000	-6.000000000000	0.707106781187
            102.236428500000	-12.562500469379	-6.000000000000	1.000000000000
            102.236428500000	-6.332138625736	-6.000000000000	0.995973631625
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000
            89.673928030621	-12.562500469379	-6.000000000000	0.707106781187
            95.904289874264	-6.332138625736	-6.000000000000	0.704259708805
            102.236428500000	-0.000000000000	-6.000000000000	0.707106781187
            89.673928030621	-0.000000000000	-6.000000000000	1.000000000000
            95.904289874264	-0.000000000000	-6.000000000000	0.995973631625
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000];
    case 'bottom_left'
        temp = [89.673928030621	-0.000000000000	-6.000000000000	1.000000000000
            95.904289874264	-0.000000000000	-6.000000000000	0.995973631625
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000
            89.673928030621	12.562500469379	-6.000000000000	0.707106781187
            95.904289874264	6.332138625736	-6.000000000000	0.704259708805
            102.236428500000	-0.000000000000	-6.000000000000	0.707106781187
            102.236428500000	12.562500469379	-6.000000000000	1.000000000000
            102.236428500000	6.332138625736	-6.000000000000	0.995973631625
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000
            114.798928969379	12.562500469379	-6.000000000000	0.707106781187
            108.568567125736	6.332138625736	-6.000000000000	0.704259708805
            102.236428500000	-0.000000000000	-6.000000000000	0.707106781187
            114.798928969379	-0.000000000000	-6.000000000000	1.000000000000
            108.568567125736	-0.000000000000	-6.000000000000	0.995973631625
            102.236428500000	-0.000000000000	-6.000000000000	1.000000000000];
    case 'top_right'
        temp = [114.798928969379 -0.000000000000 -6.879260917442 1.000000000000
            108.568567125736 -0.000000000000 -8.010000000000 0.995973631625
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000
            114.798928969379 -12.562500469379 -6.879260917442 0.707106781187
            108.568567125736 -6.332138625736 -8.010000000000 0.704259708805
            102.236428500000 -0.000000000000 -8.009999999999 0.707106781187
            102.236428500000 -12.562500469379 -6.879260917442 1.000000000000
            102.236428500000 -6.332138625736 -8.010000000000 0.995973631625
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000
            89.673928030621 -12.562500469379 -6.879260917442 0.707106781187
            95.904289874264 -6.332138625736 -8.010000000000 0.704259708805
            102.236428500000 -0.000000000000 -8.009999999999 0.707106781187
            89.673928030621 -0.000000000000 -6.879260917442 1.000000000000
            95.904289874264 -0.000000000000 -8.010000000000 0.995973631625
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000];
    case 'top_left'
        temp = [89.673928030621 -0.000000000000 -6.879260917442 1.000000000000
            95.904289874264 -0.000000000000 -8.010000000000 0.995973631625
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000
            89.673928030621 12.562500469379 -6.879260917442 0.707106781187
            95.904289874264 6.332138625736 -8.010000000000 0.704259708805
            102.236428500000 -0.000000000000 -8.009999999999 0.707106781187
            102.236428500000 12.562500469379 -6.879260917442 1.000000000000
            102.236428500000 6.332138625736 -8.010000000000 0.995973631625
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000
            114.798928969379 12.562500469379 -6.879260917442 0.707106781187
            108.568567125736 6.332138625736 -8.010000000000 0.704259708805
            102.236428500000 -0.000000000000 -8.009999999999 0.707106781187
            114.798928969379 -0.000000000000 -6.879260917442 1.000000000000
            108.568567125736 -0.000000000000 -8.010000000000 0.995973631625
            102.236428500000 -0.000000000000 -8.009999999999 1.000000000000];
end

control_point_list = PointList(temp);

nurbs = Nurbs(knot_vectors, order, control_point_list);
end

function nurbs = Solid_Cube(length, width, hight, center)
import Utility.BasicUtility.PointList
import Utility.NurbsUtility.Nurbs

knot_vectors = {[0 0 1 1] [0 0 1 1] [0 0 1 1]};
order = [1 1 1];
temp = [center(1)-length/2 center(2)-width/2 center(3)-hight/2 1.0;
        center(1)+length/2 center(2)-width/2 center(3)-hight/2 1.0;
        center(1)-length/2 center(2)+width/2 center(3)-hight/2 1.0;
        center(1)+length/2 center(2)+width/2 center(3)-hight/2 1.0;
        
        center(1)-length/2 center(2)-width/2 center(3)+hight/2 1.0;
        center(1)+length/2 center(2)-width/2 center(3)+hight/2 1.0;
        center(1)-length/2 center(2)+width/2 center(3)+hight/2 1.0;
        center(1)+length/2 center(2)+width/2 center(3)+hight/2 1.0];
    
control_point_list = PointList(temp);

nurbs = Nurbs(knot_vectors, order, control_point_list);
end

% Utilities to do knot insertion or degreee elevation
function new_nurbs = degreeElevation(nurbs, degree)
[new_points, new_knots, ~, new_order] = Nurb_DegElev(nurbs, degree);
import Utility.BasicUtility.PointList

% transform to Cartesian coordinates
for i = 1:size(new_points, 1)
    new_points(i,1:3) = new_points(i,1:3)/new_points(i,4);
end

import Utility.NurbsUtility.Nurbs

new_nurbs = Nurbs(new_knots, new_order, PointList(new_points));
end

function [Qw, Uh] = DegreeElev(p, Pw, UP, t)
% Input: p, UP, Pw, t
% Output: Uh, Qw
[N_m, N_n, N_q] = size(Pw);
n = length(UP)-p-2;
m = n+p+1;
ph = p+t;
ph2 = floor(ph/2);

%% Compute Bezier degree elevation coefficients
bezalfs(1, 1) = 1;
bezalfs(p+1, ph+1) = 1;
for i = 1:ph2
    inv = 1.0 / Bin(ph, i);
    mpi = min(p, i);
    for j = max(0, i-t):mpi
        bezalfs(j+1, i+1) = inv*Bin(p, j)*Bin(t, i-j);
    end
end
for i = (ph2+1):(ph-1)
    mpi = min(p, i);
    for j = max(0, i-t):mpi
        bezalfs(j+1, i+1) = bezalfs(p-j+1, ph-i+1);
    end
end
ua = UP(1);
for i = 0:ph
    Uh(i+1) = ua;
end

for qq = 1:N_q
    for nn = 1:N_n
        ua = UP(1);
        Local_Pw = [];
        Local_Qw=[];
        mh = ph;
        kind = ph+1;
        r = -1;
        a = p;
        b = p+1;
        cind = 1;
        for mm = 1:N_m
            Local_Pw(mm, :) = [Pw{mm, nn, qq}];
        end
        %% Initialize first Bezier seg
        [~, col_n] = size(Local_Pw);
        for col = 1:col_n
            Local_Qw(1, col) = Local_Pw(1, col);
            for i = 0:p
                bpts(i+1, col) = Local_Pw(i+1, col);
            end
        end
        % Big loop thru knot vector
        while b < m
            i = b;
            while b < m && UP(b+1) == UP(b+2)
                b = b+1;
            end
            mul = b - i + 1;
            mh = mh + mul + t;
            ub = UP(b+1);
            oldr = r;
            r = p-mul;
            
            % Insert knot u(b) r times
            if oldr > 0
                lbz = floor((oldr+2)/2);
            else
                lbz = 1;
            end
            if r > 0
                rbz = ph-floor((r+1)/2);
            else
                rbz = ph;
            end
            
            % Insert knot to get Bezier segment
            if r > 0
                numer = ub-ua;
                for k = p:-1:(mul+1)
                    alfs(k-mul) = numer / (UP(a+k+1)-ua);
                end
                for j = 1:r
                    save = r - j;
                    s = mul + j;
                    for col = 1:col_n
                        for k = p:-1:s
                            bpts(k+1,col) = alfs(k-s+1)*bpts(k+1,col) + (1.0-alfs(k-s+1))*bpts(k,col);
                        end
                        Nextbpts(save+1,col)=bpts(p+1,col);
                    end
                end
            end
            % Degree elevate Bezier
            for i = lbz:ph
                for col = 1:col_n
                    ebpts(i+1, col) = 0;
                    mpi = min(p, i);
                    for j = max(0, i-t):mpi
                        ebpts(i+1, col) = ebpts(i+1, col) + bezalfs(j+1, i+1)*bpts(j+1, col);
                    end
                end
            end
            % Must remove knot u=U(a) oldr times
            if oldr > 1
                first = kind-2;
                last = kind;
                den = ub - ua;
                bet = floor((ub-Uh(kind)) / den);
                for tr = 1:(oldr-1)
                    i = first;
                    j = last;
                    kj = j-kind+1;
                    while (j-i) > tr
                        if i < cind
                            alf = (ub-Uh(i+1)) / (ua-Uh(i+1));
                            for col = 1:col_n
                                Local_Qw(i+1, col) = alf*Local_Qw(i+1, col) + (1.0-alf)*Local_Qw(i-1+1, col);
                            end
                        end
                        if j >= lbz
                            if (j-tr) <= (kind-ph+oldr)
                                gam = (ub-Uh(j-tr+1)) / den;
                                for col = 1:col_n
                                    ebpts(kj+1, col) = gam*ebpts(kj+1, col) + (1.0-gam)*ebpts(kj+2, col);
                                end
                            else
                                for col = 1:col_n
                                    ebpts(kj+1, col) = bet*ebpts(kj+1, col) + (1.0-bet)*ebpts(kj+2, col);
                                end
                            end
                        end
                        i = i+1;
                        j = j-1;
                        kj = kj-1;
                    end
                    first = first-1;
                    last = last+1;
                end
            end
            
            % Load the knot ua
            if a ~= p
                for i = 0:(ph-oldr-1)
                    Uh(kind+1) = ua;
                    kind = kind+1;
                end
            end
            for j = lbz:rbz
                for col = 1:col_n
                    Local_Qw(cind+1, col) = ebpts(j+1, col);
                end
                cind = cind+1;
            end
            if b < m
                for col = 1:col_n
                    for j = 0:(r-1)
                        bpts(j+1, col) = Nextbpts(j+1, col);
                    end
                    for j = r:p
                        bpts(j+1, col) = Local_Pw(b-p+j+1, col);
                    end
                end
                a = b;
                b = b+1;
                ua = ub;
            else
                for i = 0:ph
                    Uh(kind+i+1) = ub;
                end
            end
        end % End big while loop
        
        [rowQ, ~] = size(Local_Qw);
        for i = 1:rowQ
            Qw{i, nn, qq} = Local_Qw(i, :);
        end
    end
end
    function b = Bin(p,i)
        b = factorial(p) / (factorial(i)*factorial(p-i));
    end
end

function [new_points, new_knots, new_basis_number, new_order] = Nurb_DegElev(nurbs, t)
% Input: p(matrix), knots(cell), point(matrix), t(matrix), basis_num(matrix)
% Output: new_points, new_knots, new_basis_number, new_order
import Utility.BasicUtility.TensorProduct

knots = nurbs.knot_vectors_;
p = nurbs.order_;
point = nurbs.control_points_(:,:);
basis_n = nurbs.basis_number_;

% transform to homogeneous coordinates
for i = 1:size(point,1)
    point(i,1:3) = point(i,1:3)*point(i,4);
end

while length(basis_n) < 3
    basis_n(1, end+1) = 1;
end

%Rearrange point list to cubic
TD = TensorProduct({basis_n(1) basis_n(2) basis_n(3)});
point_L = basis_n(1) * basis_n(2) * basis_n(3);
for k = 1:point_L
    pos = TD.to_local_index(k);
    Pw(pos{1},pos{2},pos{3}) = {point(k,:)};
end

Dim = length(knots);
switch Dim
    case 1
        UP = knots{1};
        t_u = t(1);
        %insert u-direction
        if t_u == 0
            Uh = UP;
        else
            [Qw, Uh] = DegreeElev(p(1), Pw, UP, t_u);
            Pw = Qw;
        end
        new_knots = {Uh};
        
    case 2
        UP = knots{1}; VP = knots{2};
        t_u = t(1); t_v = t(2);
        %elevate u-direction
        if t_u == 0
            Uh = UP;
        else
            [Qw, Uh] = DegreeElev(p(1), Pw, UP, t_u);
            Pw = Qw;
        end
        %elevate v-direction
        if t_v == 0
            Vh = VP;
        else
            tempPw = permute(Pw, [2 1 3]);
            [Qw, Vh] = DegreeElev(p(2), tempPw, VP, t_v);
            Pw = permute(Qw,[2 1 3]);
        end
        new_knots = {Uh, Vh};
        
    case 3
        UP = knots{1}; VP = knots{2}; WP = knots{3};
        t_u = t(1); t_v = t(2); t_w = t(3);
        %insert u-direction
        if t_u == 0
            Uh = UP;
        else
            [Qw, Uh] = DegreeElev(p(1), Pw, UP, t_u);
            Pw = Qw;
        end
        %elevate v-direction
        if t_v == 0
            Vh = VP;
        else
            tempPw = permute(Pw, [2 1 3]);
            [Qw, Vh] = DegreeElev(p(2), tempPw, VP, t_v);
            Pw = permute(Qw,[2 1 3]);
        end
        %elevate w-direction
        if t_w == 0
            Wh = WP;
        else
            tempPw = permute(Pw,[3 2 1]);
            [Qw, Wh] = DegreeElev(p(3), tempPw, WP, t_w);%
            Pw = permute(Qw,[3 2 1]);
        end
        new_knots = {Uh, Vh, Wh};
end

%Rearrange cubic control point to point-list type
[N_i, N_j, N_k] = size(Pw);
TD_3 = TensorProduct({N_i N_j N_k});
for k = 1:N_k
    for j = 1:N_j
        for i = 1:N_i
            global_index = TD_3.to_global_index({i j k});
            new_points(global_index, :) = Pw{i, j, k};
        end
    end
end
new_basis_number = size(Pw);
new_order = nurbs.order_ + t;
end

function new_nurbs = knotInsertion(nurbs, knots)
[new_points, new_knots, ~] = Nurb_KnotIns(nurbs, knots);
import Utility.BasicUtility.PointList

% transform to Cartesian coordinates
for i = 1:size(new_points, 1)
    new_points(i,1:3) = new_points(i,1:3)/new_points(i,4);
end

import Utility.NurbsUtility.Nurbs

new_nurbs = Nurbs(new_knots, nurbs.order_, PointList(new_points));

end

function [Qw, UQ] = KnotIns(p, Pw, UP, u, r)
% Input: p, UP, Pw, u, r,
% Output: UQ, Qw
import BasisFunction.IGA.NurbsBasisFunction

[N_m, N_n, N_q] = size(Pw);
np = length(UP)-p-2;
k = NurbsBasisFunction.FindSpan(np, p, u, UP);
k = k+1;
mp = np+p+1;
nq = np+r;

UQ=zeros(1,length(UP)+r);

for i = 1:k
    UQ(i) = UP(i);
end
for i = 1:r
    UQ(k+i) = u;
end
for i = (k+1):(mp+1)
    UQ(i+r) = UP(i);
end

for j = 1:r
    L = k-p+j;
    for i = 0:(p-j)
        alpha(i+1, j) = (u-UP(L+i))/(UP(i+1+k)-UP(L+i));
    end
end

for q = 1:N_q
    for n = 1:N_n
        Local_Pw = [];
        for m = 1:N_m
            Local_Pw(m, :) = [Pw{m, n, q}];
        end
        [row_n, col_n] = size(Local_Pw);
        for col = 1:col_n
            for i = 0:(k-p-1)
                Local_Qw(i+1, col) = Local_Pw(i+1, col);
            end
            for i = k:(np+1)
                Local_Qw(i+r, col) = Local_Pw(i, col);
            end
            for i = 0:p
                Rw(i+1) = Local_Pw(k-p+i, col);
            end
            for j = 1:r
                L = k-p+j;
                for i = 0:(p-j)
                    Rw(i+1) = alpha(i+1,j) * Rw(i+2) + (1.0-alpha(i+1,j)) * Rw(i+1);
                    Local_Qw(L, col) = Rw(1);
                    Local_Qw(k+r-j, col) = Rw(p-j+1);
                end
                for i = (L+1):k
                    Local_Qw(i, col) = Rw(i-L+1);
                end
            end
        end
        [rowQ,colQ] = size(Local_Qw);
        for i = 1:rowQ
            Qw{i, n, q} = Local_Qw(i, :);
        end
    end
end
end

function [new_points, new_knots, new_basis_number] = Nurb_KnotIns(nurbs, u)
% Input: p(matrix), knots(cell), point(matrix), u(cell), basis_num(cell)
% Output: new_points, new_knots
import Utility.BasicUtility.TensorProduct

knots = nurbs.knot_vectors_;
p = nurbs.order_;
point = nurbs.control_points_(:,:);

% transform to homogeneous coordinates
for i = 1:size(point,1)
    point(i,1:3) = point(i,1:3)*point(i,4);
end

basis_n = nurbs.basis_number_;

while length(basis_n) < 3
    basis_n(1,end+1) = 1;
end

%Rearrange point list to cubic
TD = TensorProduct({basis_n(1) basis_n(2) basis_n(3)});
point_L = basis_n(1) * basis_n(2) * basis_n(3);
for k = 1:point_L
    pos = TD.to_local_index(k);
    Pw(pos{1},pos{2},pos{3}) = {point(k,:)};
end

Dim = length(knots);
switch Dim
    case 1
        UP = knots{1};
        u_u = u{1};
        %insert u-durection
        if isempty(u_u)
            UQ = UP;
        else
            uni = unique(u_u);
            n = histc(u_u,uni);
            Qw = Pw;
            U = UP;
            for i = 1:length(uni)
                local_u = uni(i);
                r = n(i);
                [Qw, U] = KnotIns(p(1), Qw, U, local_u, r); %
            end
            UQ = U;
            Pw = Qw;
        end
        new_knots = {UQ};
        
    case 2
        UP = knots{1}; VP = knots{2};
        u_u = u{1}; v_u = u{2};
        %insert u-durection
        if isempty(u_u)
            UQ = UP;
        else
            uni = unique(u_u);
            n = histc(u_u,uni);
            Qw = Pw;
            U = UP;
            for i = 1:length(uni)
                local_u = uni(i);
                r = n(i);
                [Qw, U] = KnotIns(p(1), Qw, U, local_u, r); %
            end
            UQ = U;
            Pw = Qw;
        end
        %insert v-durection
        if isempty(v_u)
            VQ = VP;
        else
            uni = unique(v_u);
            n = histc(v_u,uni);
            Qw = permute(Pw, [2 1 3]);
            U = VP;
            for i = 1:length(uni)
                local_u = uni(i);
                r = n(i);
                [Qw, U] = KnotIns(p(2), Qw, U, local_u, r); %
            end
            VQ = U;
            Pw = permute(Qw,[2 1 3]);
        end
        new_knots = {UQ, VQ};
        
    case 3
        UP = knots{1}; VP = knots{2}; WP = knots{3};
        u_u = u{1}; v_u = u{2}; w_u = u{3};
        %insert u-durection
        if isempty(u_u)
            UQ = UP;
        else
            uni = unique(u_u);
            n = histc(u_u,uni);
            Qw = Pw;
            U = UP;
            for i = 1:length(uni)
                local_u = uni(i);
                r = n(i);
                [Qw, U] = KnotIns(p(1), Qw, U, local_u, r); %
            end
            UQ = U;
            Pw = Qw;
        end
        
        %insert v-durection
        if isempty(v_u)
            VQ = VP;
        else
            uni = unique(v_u);
            n = histc(v_u,uni);
            Qw = permute(Pw,[2 1 3]);
            U = VP;
            for i = 1:length(uni)
                local_u = uni(i);
                r = n(i);
                [Qw, U] = KnotIns(p(2), Qw, U, local_u, r); %
            end
            VQ = U;
            Pw = permute(Qw,[2 1 3]);
        end
        
        %insert w-durection
        if isempty(w_u)
            WQ = WP;
        else
            uni = unique(w_u);
            n = histc(w_u,uni);
            Qw = permute(Pw,[3 2 1]);
            U = WP;
            for i = 1:length(uni)
                local_u = uni(i);
                r = n(i);
                [Qw, U] = KnotIns(p(3), Qw, U, local_u, r);%
            end
            WQ = U;
            Pw = permute(Qw,[3 2 1]);
        end
        new_knots = {UQ, VQ, WQ};
end

%Rearrange cubic control point to point-list type
[N_i, N_j, N_k] = size(Pw);
TD_3 = TensorProduct({N_i N_j N_k});
for k = 1:N_k
    for j = 1:N_j
        for i = 1:N_i
            global_index = TD_3.to_global_index({i j k});
            new_points(global_index, :) = Pw{i, j, k};
        end
    end
end

new_basis_number = size(Pw);
end



