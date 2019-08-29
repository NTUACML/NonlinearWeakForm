function DemoIGA_Uniaxial_Tension_Beam
clc; clear; close all;

%% Include package
import Utility.BasicUtility.*
import Utility.NurbsUtility.*
import Geometry.*
import Domain.*
import Operation.*
import Interpolation.IGA.Interpolation

%% Geometry data input
xml_path = './ArtPDE_IGA_Rectangle.art_geometry';
geo = GeometryBuilder.create('IGA', 'XML', xml_path);
nurbs_topology = geo.topology_data_{1};

%% Domain create
iga_domain = DomainBuilder.create('IGA');

%% Basis create
nurbs_basis = iga_domain.generateBasis(nurbs_topology);

%% Nurbs tools create & plot nurbs
nurbs_tool = NurbsTools(nurbs_basis);

figure; hold on; grid on; axis equal;
nurbs_tool.plotNurbs();
nurbs_tool.plotControlMesh();

hold off;

%% Variable define   
var_u = iga_domain.generateVariable('displacement', nurbs_basis,...
                                    VariableType.Vector, 2);      
%% Test variable define
test_u = iga_domain.generateTestVariable(var_u, nurbs_basis);

%% Operation define (By User)
operation1 = Operation();
operation1.setOperator('delta_epsilon_dot_sigma');

% operation2 = Operation();
% operation2.setOperator('elasticity_nitsche_dirichlet_lhs_term');

operation3 = Operation();
operation3.setOperator('elasticity_applied_traction');

%% Expression acquired
import MaterialBank.*
import Utility.BasicUtility.ElasticMaterialType
E = 2e8;
nu = 0.3;
constitutive_law = SaintVenantKirchhoff({nu, E, ElasticMaterialType.PlaneStress});

exp1 = operation1.getExpression('IGA', {test_u, var_u, constitutive_law});

% beta = 100;
% exp2 = operation2.getExpression('IGA', {test_u, var_u, constitutive_law, beta});

exp3 = operation3.getExpression('IGA', {test_u, @traction});

%% Integral variation equations
% Domain integral
import Differential.IGA.*

% Domain integral
domain_patch = nurbs_topology.getDomainPatch();
dOmega = Differential(nurbs_basis, domain_patch);

iga_domain.integrate(exp1, dOmega, {'Default', 2});

bdr_patch = nurbs_topology.getBoundayPatch('xi_1');
dGamma = Differential(nurbs_basis, bdr_patch);

iga_domain.integrate(exp3, dGamma, {'Default', 2});

%% Constraint (Acquire prescribed D.O.F.)
bdr_patch = nurbs_topology.getBoundayPatch('eta_0');
iga_domain.generateConstraint(bdr_patch, var_u, {2, @zero_fcn}, {'collocation', nurbs_basis});

bdr_patch = nurbs_topology.getBoundayPatch('xi_0');
iga_domain.generateConstraint(bdr_patch, var_u, {1, @zero_fcn}, {'collocation', nurbs_basis});
%% Solve domain equation system
import Solver.*
import Utility.BasicUtility.SolverType
solver = Solver();
solver.generate(iga_domain, SolverType.Standard);
solver.solve();

%% Data Interpolation
t_interpo = Interpolation(var_u);
[x, data, element] = t_interpo.DomainDataSampling(11, 'gradient');

%% Show result (Post-Processes)
% plot x-displacement
fv.vertices = [x(:,1:2), data.value{1}];
fv.faces = element;
fv.facevertexcdata = data.value{1};

figure; hold on; grid on; view([0 90]);
patch(fv,'CDataMapping','scaled','EdgeColor',[.7 .7 .7],'FaceColor','interp','FaceAlpha',1);
title('IGA Uniaxial Tension Beam x-displacement')
colorbar
hold off;

% plot y-displacement
fv.vertices = [x(:,1:2), data.value{2}];
fv.faces = element;
fv.facevertexcdata = data.value{2};

figure; hold on; grid on; view([0 90]);
patch(fv,'CDataMapping','scaled','EdgeColor',[.7 .7 .7],'FaceColor','interp','FaceAlpha',1);
title('IGA Uniaxial Tension Beam y-displacement')
colorbar
hold off;

% plot deformed mesh
import Utility.Resources.quadplot
figure; hold on; grid on; axis equal;
scale_factor = 1e5;
quadplot(element, x(:,1)+scale_factor*data.value{1}, x(:,2)+scale_factor*data.value{2});
title('IGA Cantilever Beam deformed mesh')
hold off;

temp = E/(1-nu*nu);
% plot sigma_x
s_x = (data.gradient{1,1} + nu*data.gradient{2,2})*temp;

fv.vertices = [x(:,1:2), s_x];
fv.faces = element;
fv.facevertexcdata = s_x;

figure; hold on; grid on; view([20 30]);
patch(fv,'CDataMapping','scaled','EdgeColor',[.7 .7 .7],'FaceColor','interp','FaceAlpha',1);
title('IGA Cantilever Beam \sigma_x')
view([20 30]);

ana = sigma_x(x(:,1:2));
plot3(x(:,1), x(:,2), ana, 'r.');
rms_ = rms(s_x-ana);
disp(rms_);
hold off;

% plot sigma_y
s_y = (nu*data.gradient{1,1} + data.gradient{2,2})*temp;

fv.vertices = [x(:,1:2), s_y];
fv.faces = element;
fv.facevertexcdata = s_y;

figure; hold on; grid on; view([20 30]);
patch(fv,'CDataMapping','scaled','EdgeColor',[.7 .7 .7],'FaceColor','interp','FaceAlpha',1);
title('IGA Cantilever Beam \sigma_y')
view([20 30]);

ana = sigma_y(x(:,1:2));
plot3(x(:,1), x(:,2), ana, 'r.');
rms_ = rms(s_y-ana);
disp(rms_);
hold off;

% plot sigma_xy
s_xy = (data.gradient{1,2} + data.gradient{2,1})*temp*0.5*(1-nu);

fv.vertices = [x(:,1:2), s_xy];
fv.faces = element;
fv.facevertexcdata = s_xy;

figure; hold on; grid on; 
patch(fv,'CDataMapping','scaled','EdgeColor',[.7 .7 .7],'FaceColor','interp','FaceAlpha',1);
title('IGA Cantilever Beam \sigma_{xy}')
view([20 30]);

ana = sigma_xy(x(:,1:2));
plot3(x(:,1), x(:,2), ana, 'r.');
rms_ = rms(s_xy-ana);
disp(rms_);
hold off;

%% Show result
disp(var_u);
end

function [t_x, t_y] = traction(x, y)
    t_x = 1e3;
    t_y = 0;
end

function val = zero_fcn(x, y)
    val = 0;
end
