function DemoIGA_LinearElasticity3D
clc; clear; close all;
%% Include package
import Utility.BasicUtility.*
import Utility.NurbsUtility.*
import Geometry.*
import Domain.*
import Operation.*
import Interpolation.IGA.Interpolation

%% Geometry data input
xml_path = './ArtPDE_IGA_Unit_Cube.art_geometry';
geo = GeometryBuilder.create('IGA', 'XML', xml_path);
nurbs_topology = geo.topology_data_{1};

%% Domain create
iga_domain = DomainBuilder.create('IGA');

%% Basis create
nurbs_basis = iga_domain.generateBasis(nurbs_topology);

%% Nurbs tools create & plot nurbs
nurbs_tool = NurbsTools(nurbs_basis);

figure; hold on; grid on; view([120 30]); %axis equal;
nurbs_tool.plotNurbs();
xlabel('x'); ylabel('y'); zlabel('z'); 

% nurbs_tool.plotControlMesh();
control_point = nurbs_topology.domain_patch_data_.nurbs_data_.control_points_(:,1:3);
for i = 1:size(control_point,1)
    text(control_point(i,1), control_point(i,2), control_point(i,3), num2str(i), 'FontSize',14);
end

% hold off;

%% Variable define   
u = iga_domain.generateVariable('displacement', nurbs_basis,...
                                    VariableType.Vector, 3);                               
%% Test variable define
v = iga_domain.generateTestVariable(u, nurbs_basis);

%% Operation define (By User)
operation1 = Operation();
operation1.setOperator('ElasticityBilinearExpression3D');

operation2 = Operation();
operation2.setOperator('ElasticityTractionExpression3D');

%% Expression acquired
import MaterialBank.*
import Utility.BasicUtility.ElasticMaterialType

E = 1e7;
nu = 0.3;
constitutive_law = SaintVenantKirchhoff({nu, E, ElasticMaterialType.Solid});

exp1 = operation1.getExpression('IGA', {v, u, constitutive_law});

exp2 = operation2.getExpression('IGA', {v, @traction});
%% Integral variation equations
% Domain integral
import Differential.IGA.*

% Domain integral
domain_patch = nurbs_topology.getDomainPatch();
dOmega = Differential(nurbs_basis, domain_patch);

iga_domain.integrate(exp1, dOmega, {'Default', 2});

% Traction integral
bdr_patch = nurbs_topology.getBoundayPatch('xi_1');
dGamma = Differential(nurbs_basis, bdr_patch);

iga_domain.integrate(exp2, dGamma, {'Default', 2});

%% Constraint (Acquire prescribed D.O.F.)
bdr_patch = nurbs_topology.getBoundayPatch('xi_0');
iga_domain.generateConstraint(bdr_patch, u, {1, @()0});

bdr_patch = nurbs_topology.getBoundayPatch('eta_0');
iga_domain.generateConstraint(bdr_patch, u, {2, @()0});

bdr_patch = nurbs_topology.getBoundayPatch('zeta_0');
iga_domain.generateConstraint(bdr_patch, u, {3, @()0});

% temp = iga_domain.constraint_(3);
% temp.constraint_var_id_ = [1 3];
%% Solve domain equation system
import Solver.*
import Utility.BasicUtility.SolverType
solver = Solver();

solver.generate(iga_domain, SolverType.Standard);
solver.solve();

%% Show result
disp(u);

%% Show result (Post-Processes)
% plot deformed nurbs
scaleFactor = 10000;
nurbs_topology.domain_patch_data_.nurbs_data_.control_points_(:,1:3) = ...
    nurbs_topology.domain_patch_data_.nurbs_data_.control_points_(:,1:3) + scaleFactor*u.getVarData();

nurbs_tool.plotNurbs();
hold off;


end

function [t_x, t_y, t_z] = traction(x, y, z)
    t_x = 1/(0.1*0.2);
    t_y = 0;
    t_z = 0;
end

