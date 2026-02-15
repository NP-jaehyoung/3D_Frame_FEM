%................................................................
% MATLAB codes for Finite Element Analysis
% problem13.m
% antonio ferreira 2008
% clear memory
clear all
clc
% E; modulus of elasticity [unit:kPa]
% I: second moments of area
% J: polar moment of inertia
% G: shear modulus
% L: length of bar
% unit : force : kN / Disp : m 
E=210e6; A=0.02; 
Iy=10e-5;   Iz=20e-5; J=5e-5; G=84e6;
% generation of coordinates and connectivities
nodeCoordinates=[0 0 0; %1st node
    0 0 4; %2nd node
    4 0 4; %3rd node
    4 0 0; %4th node

    0 5 0; %5th node
    0 5 4; %6th node
    4 5 4; %7th node
    4 5 0;]; %8th node

xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
zz=nodeCoordinates(:,3);
% elementNodes=setEleNodes(nodeCoordinates);
elementNodes=[1 5; % Column
    2 6; % Column
    3 7; % Column
    4 8; % Column
    5 6; % Girder
    6 7; % Girder
    7 8; % Girder
    8 5;]; % Girder 
numberNodes=size(nodeCoordinates,1);
numberElements=size(elementNodes,1);
% for structure:
% displacements: displacement vector
% force : force vector
% stiffness: stiffness matrix
% GDof: global number of degrees of freedom
GDof=6*numberNodes;
U=zeros(GDof,1);
force=zeros(GDof,1);
stiffness=zeros(GDof);
%force vector
force(37)=-1500;
% calculation of the system stiffness matrix
% and force vector
% stiffness matrix
[stiffness]=...
    formStiffness3Dframe(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,E,A,Iz,Iy,G,J);
% boundary conditions and solution
prescribedDof=[1:24];
% solution
displacements=solution(GDof,prescribedDof,stiffness,force);
% displacements
disp("Displacements")
jj=1:GDof; format long
f=[jj; displacements'];
fprintf("node U\n")
fprintf("%3d %12.8f\n",f)
%drawing mesh and deformed shape
U=displacements;
dispScale = 5;
clf
figure(2)
hold on
drawingMesh(nodeCoordinates+dispScale*[U(1:6:6*numberNodes) U(2:6:6*numberNodes) U(3:6:6*numberNodes)], elementNodes,"L2","k.-");
hold on
drawingMesh(nodeCoordinates,elementNodes,"L2","k--");
title("Scaled Disp")
xlabel("X-Coord")
ylabel("Y-Coord")
zlabel("Z-Coord")
legend("Deform","UnDeform")
function  [stiffness]=formStiffness3Dframe(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,E,A,Iz,Iy,G,J);
stiffness=zeros(GDof);
% computation of the system stiffness matrix
for e=1:numberElements
    % elementDof: element degrees of freedom (Dof)
    indice=elementNodes(e,:)   ;
    elementDof=[6*indice(1)-5 6*indice(1)-4 6*indice(1)-3 ...
        6*indice(1)-2 6*indice(1)-1 6*indice(1)...
        6*indice(2)-5 6*indice(2)-4 6*indice(2)-3 ...
        6*indice(2)-2 6*indice(2)-1 6*indice(2)] ;
    x1=nodeCoordinates(indice(1),1);
    y1=nodeCoordinates(indice(1),2);
    z1=nodeCoordinates(indice(1),3);
    x2=nodeCoordinates(indice(2),1);
    y2=nodeCoordinates(indice(2),2);
    z2=nodeCoordinates(indice(2),3);
    L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +...
        (z2-z1)*(z2-z1));
    k1 = E*A/L;
    k2 = 12*E*Iz/(L*L*L);
    k3 = 6*E*Iz/(L*L);
    k4 = 4*E*Iz/L;
    k5 = 2*E*Iz/L;
    k6 = 12*E*Iy/(L*L*L);
    k7 = 6*E*Iy/(L*L);
    k8 = 4*E*Iy/L;
    k9 = 2*E*Iy/L;
    k10 = G*J/L;
    a=[k1 0 0; 0 k2 0; 0 0 k6];
    b=[ 0 0 0;0 0 k3; 0 -k7 0];
    c=[k10 0 0;0 k8 0; 0 0 k4];
    d=[-k10 0 0;0 k9 0;0 0 k5];
    k = [a b -a b;b' c b d; (-a)' b' a -b;b' d' (-b)' c];
    if x1 == x2 && y1 == y2
        if z2 > z1
            Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
        else
            Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
        end
    else


        CXx = (x2-x1)/L;
        CYx = (y2-y1)/L;
        CZx = (z2-z1)/L;
        D = sqrt(CXx*CXx + CYx*CYx);
        CXy = -CYx/D;
        CYy = CXx/D;
        CZy = 0;
        CXz = -CXx*CZx/D;
        CYz = -CYx*CZx/D;
        CZz = D;
        Lambda = [CXx CYx CZx ;CXy CYy CZy ;CXz CYz CZz];
    end
    R = [Lambda zeros(3,9); zeros(3) Lambda zeros(3,6);
        zeros(3,6) Lambda zeros(3);zeros(3,9) Lambda];
    stiffness(elementDof,elementDof)=...
        stiffness(elementDof,elementDof)+R'*k*R;
end
end

function displacements=solution(GDof,prescribedDof,stiffness,force)
% function to find solution in terms of global displacements
activeDof=setdiff([1:GDof]', [prescribedDof]);
U=stiffness(activeDof,activeDof)\force(activeDof);
displacements=zeros(GDof,1);
displacements(activeDof)=U;
end