%................................................................
% MATLAB codes for Finite Element Analysis
% problem13.m
% antonio ferreira 2008
% clear memory
clear all
clc
% E; modulus of elasticity [unit:kPa]
% I: second moments of area
% J: polar moment of inertia /// 24.12.18. AN ; tortional resistance??
% G: shear modulus
% L: length of bar
% unit : force : kN / Disp : m
E=2e7; width = 1; depth = 1; A=width*depth;
nu=0.167 
Iy=1/12;  
Iz=1/12;
J=2.25*1^4; 
G=E/(2*(1+nu));  % G=E/(2(1+ν))  ν:포아송비
% generation of coordinates and connectivities
nodeCoordinates=[0 0 0; %1st node
    0 0 4; %2nd node
    4 0 4; %3rd node
    4 0 0; %4th node

    0 5 0; %5th node
    0 5 4; %6th node
    4 5 4; %7th node
    4 5 0;]; %8th node
%%% stories you want %%%
numFloors = 2;
for i=3:numFloors
    nodeCoordinates(4*(i-1)+1,:) = [0 5*(i-1) 0];
    nodeCoordinates(4*(i-1)+2,:) = [0 5*(i-1) 4];
    nodeCoordinates(4*(i-1)+3,:) = [4 5*(i-1) 4];
    nodeCoordinates(4*(i-1)+4,:) = [4 5*(i-1) 0];
end
%%
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
for i=3:numFloors
    elementNodes(8*(i-2)+1,:) = [4*(i-2)+1 4*(i-1)+1]; % Column
    elementNodes(8*(i-2)+2,:) = [4*(i-2)+2 4*(i-1)+2]; % Column
    elementNodes(8*(i-2)+3,:) = [4*(i-2)+3 4*(i-1)+3]; % Column
    elementNodes(8*(i-2)+4,:) = [4*(i-2)+4 4*(i-1)+4]; % Column
    elementNodes(8*(i-2)+5,:) = [4*(i-1)+1 4*(i-1)+2]; % Girder
    elementNodes(8*(i-2)+6,:) = [4*(i-1)+2 4*(i-1)+3]; % Girder
    elementNodes(8*(i-2)+7,:) = [4*(i-1)+3 4*(i-1)+4]; % Girder
    elementNodes(8*(i-2)+8,:) = [4*(i-1)+4 4*(i-1)+1]; % Girder
end
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
force(37)=-3000; % 37-th DOF는 7-th Node의Dx
% calculation of the system stiffness matrix
% and force vector
% stiffness matrix
[stiffness]=...
    formStiffness3Dframe(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,E,A,Iz,Iy,G,J);
% boundary conditions and solution
temp1=stiffness;
%stiffness = floorConstraint(GDof,stiffness,numFloors,nodeCoordinates);
%stiffness-temp1
prescribedDof=[1:12];
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
dispScale = 10;
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

function [stiffness] = formStiffness2Dframe(GDof, numberElements, elementNodes, numberNodes, nodeCoordinates, E, A, I)
% 이 함수는 2차원 프레임 요소 (u_x, u_y, r_z 자유도) 기반 
% 전체 구조물의 전역 강성행렬을 구성합니다.

% 초기화
stiffness = zeros(GDof);

% 각 요소에 대해 반복
for e = 1:numberElements
    % 요소를 구성하는 두 노드
    nodes = elementNodes(e,:);
    n1 = nodes(1);
    n2 = nodes(2);
    
    % 노드 좌표 추출
    x1 = nodeCoordinates(n1,1); y1 = nodeCoordinates(n1,2);
    x2 = nodeCoordinates(n2,1); y2 = nodeCoordinates(n2,2);
    
    % 요소 길이와 회전각도
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    cosT = (x2 - x1)/L;
    sinT = (y2 - y1)/L;
    
    % 요소 DOF 번호 설정
    % 각 노드는 u_x, u_y, r_z 순서로 3개 DOF
    % node n에 대한 DOF: [3*n-2, 3*n-1, 3*n]
    elementDof = [3*n1-2 3*n1-1 3*n1 3*n2-2 3*n2-1 3*n2];
    
    % 국소 강성계수
    EA = E*A;
    EI = E*I;
    
    % 국소 강성행렬(K_local)
    kLocal = [ EA/L      0            0      -EA/L      0            0;
               0         12*EI/L^3    6*EI/L^2 0       -12*EI/L^3   6*EI/L^2;
               0         6*EI/L^2     4*EI/L   0       -6*EI/L^2    2*EI/L;
              -EA/L      0            0       EA/L      0            0;
               0        -12*EI/L^3   -6*EI/L^2 0        12*EI/L^3   -6*EI/L^2;
               0         6*EI/L^2     2*EI/L   0       -6*EI/L^2    4*EI/L ];
           
    % 회전행렬 T (3x3)
    T = [cosT -sinT 0;
         sinT  cosT 0;
         0     0    1];
    
    % 블록 대각형 형태의 R (6x6)
    R = [T zeros(3,3);
         zeros(3,3) T];
    
    % 전역 좌표계로 변환
    kGlobal = R' * kLocal * R;
    
    % 전역 강성행렬에 조립
    stiffness(elementDof, elementDof) = stiffness(elementDof, elementDof) + kGlobal;
end

end


function displacements=solution(GDof,prescribedDof,stiffness,force)
% function to find solution in terms of global displacements
activeDof=setdiff([1:GDof]', [prescribedDof]);
U=stiffness(activeDof,activeDof)\force(activeDof);
displacements=zeros(GDof,1);
displacements(activeDof)=U;
end

function stiff=floorConstraint(GDof,stiffness,numFloors,nodeCoordinates)
for floor = 2:numFloors
    % Find nodes at the current floor height
    floorHeight = (floor - 1) * 5; % Example: 5m between floors
    floorNodes = find(nodeCoordinates(:,2) == floorHeight);
        % Select the first node as the master node for simplicity
    masterNode = floorNodes(1)
    
    % Link DOFs of slave nodes to the master node
    for i = 2:length(floorNodes)
        slaveNode = floorNodes(i)
        
        % Link horizontal displacements (X and Z) and rotation about Y
        % Assuming the DOFs are ordered [u1, v1, w1, theta_x1, theta_y1, theta_z1, u2, ..., theta_zn]
        % where u, v, w are displacements in x, y, z and theta_x, theta_y, theta_z are rotations
        
        % Displacement in X
        stiffness(6*masterNode-5, 6*slaveNode-5) = stiffness(6*masterNode-5, 6*slaveNode-5) + 1e5; % X-disp
        stiffness(6*slaveNode-5, 6*masterNode-5) = stiffness(6*slaveNode-5, 6*masterNode-5) + 1e5; % Symetric entty
        
        % Displacement in Z
        stiffness(6*masterNode-3, 6*slaveNode-3) = stiffness(6*masterNode-3, 6*slaveNode-3) + 1e5; % Z-disp
        stiffness(6*slaveNode-3, 6*masterNode-3) = stiffness(6*slaveNode-3, 6*masterNode-3) + 1e5; % Symetric entty

        % Rotation about Y
        stiffness(6*masterNode-2, 6*slaveNode-2) = stiffness(6*masterNode-2, 6*slaveNode-2) + 1e5; % Y-Rot
        stiffness(6*slaveNode-2, 6*masterNode-2) = stiffness(6*slaveNode-2, 6*masterNode-2) + 1e5; % Symetric entty
    end

end
stiff=stiffness;
end





