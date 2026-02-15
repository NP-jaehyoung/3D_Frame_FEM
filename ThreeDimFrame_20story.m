%................................................................
% MATLAB codes for Finite Element Analysis
% 20-story 3D frame example
% clear memory
clear all
clc

%% model settings
numFloors = 20;            % number of floor levels
floorHeight = 4.0;         % story height [m]
spanX = 4.0;               % x direction span [m]
spanZ = 4.0;               % z direction span [m]
E = 2e7;                   % modulus of elasticity [kPa]
rho = 7850;                % density
width = 1;                 
depth = 1;                 
A = width*depth;           % area [m^2]
Iy = 1/12;
Iz = 1/12;
nu = 0.167;
J = 2.25*1^4;             
G = E/(2*(1+nu));
Ptop = -3000;              % single lateral load at top story [kN]

%% node coordinates: first two floors
nodeCoordinates = [...
    0       0          0;  % 1
    0       0        spanZ; % 2
    spanX  0        spanZ; % 3
    spanX  0          0;  % 4
    0       floorHeight  0; % 5
    0       floorHeight spanZ; % 6
    spanX  floorHeight spanZ; % 7
    spanX  floorHeight  0; % 8
];
for i = 3:numFloors
    baseY = floorHeight*(i-1);
    nodeCoordinates(4*(i-1)+1,:) = [0 baseY 0];
    nodeCoordinates(4*(i-1)+2,:) = [0 baseY spanZ];
    nodeCoordinates(4*(i-1)+3,:) = [spanX baseY spanZ];
    nodeCoordinates(4*(i-1)+4,:) = [spanX baseY 0];
end

%% element connectivity: columns + girders each floor
elementNodes = [1 5; 2 6; 3 7; 4 8; 5 6; 6 7; 7 8; 8 5];
for i = 3:numFloors
    base = 4*(i-2);
    elementNodes(8*(i-2)+1,:) = [base+1 base+5];
    elementNodes(8*(i-2)+2,:) = [base+2 base+6];
    elementNodes(8*(i-2)+3,:) = [base+3 base+7];
    elementNodes(8*(i-2)+4,:) = [base+4 base+8];
    elementNodes(8*(i-2)+5,:) = [base+5 base+6];
    elementNodes(8*(i-2)+6,:) = [base+6 base+7];
    elementNodes(8*(i-2)+7,:) = [base+7 base+8];
    elementNodes(8*(i-2)+8,:) = [base+8 base+5];
end

numberNodes = size(nodeCoordinates,1);
numberElements = size(elementNodes,1);
GDof = 6*numberNodes;

%% assemble K, M
force = zeros(GDof,1);
% top story one-point load at UX of node 3 of top floor
topNode = 4*(numFloors-1) + 3;   % 1st corner node on top floor
topLoadDof = 6*topNode - 5;      % UX at that node
force(topLoadDof) = Ptop;

stiffness = formStiffness3Dframe(GDof, numberElements, elementNodes,...
    numberNodes, nodeCoordinates, E, A, Iz, Iy, G, J);
mass = formMass3Dframe(GDof, numberElements, elementNodes,...
    numberNodes, nodeCoordinates, rho, A, Iz, Iy);
stiffnessAll = stiffness;
massAll = mass;
forceAll = force;

%% rigid diaphragm constraints: tie all nodes by floor in UX, UZ, RY (dx, dz, ry)
[stiffness, mass, force, T, dofToReduced] = ...
    applyRigidDiaphragmConstraints(stiffnessAll, massAll, forceAll, ...
    nodeCoordinates, 2, numFloors);

%% boundary conditions: fix 1st floor nodes (nodes 1~4)
prescribedDof = 1:24;
prescribedDofReduced = unique(dofToReduced(prescribedDof));

%% static solution and output
displacementsReduced = solution(size(stiffness,1), prescribedDofReduced, ...
    stiffness, force);
displacements = T * displacementsReduced;
reactions = stiffnessAll * displacements - forceAll;
reactionSupport = reactions(prescribedDof);

fprintf("20-floor 3D frame example\n");
fprintf("Top node: %d, top load DOF: %d, load: %.4g\n", topNode, topLoadDof, Ptop);
fprintf("Top node displacement (UX,UY,UZ) [m]\n");
fprintf("UX = %.6e\n", displacements(topLoadDof));
fprintf("UY = %.6e\n", displacements(topLoadDof+1));
fprintf("UZ = %.6e\n\n", displacements(topLoadDof+2));

fprintf("Support reactions (first 24 DOF)\n");
for ii = 1:numel(prescribedDof)
    fprintf("DOF %d : %+13.6e\n", prescribedDof(ii), reactionSupport(ii));
end

%% modal analysis
maxModesToShow = 10;
[modesReduced,eigenvalues] = eigenvalue(size(stiffness,1), ...
    prescribedDofReduced, stiffness, mass, maxModesToShow);
modes = T * modesReduced;
omega = sqrt(abs(eigenvalues));           % rad/s
freqHz = omega/(2*pi);                   % Hz
[freqHzSorted, sortIdx] = sort(freqHz,"ascend");
fprintf("\nNatural frequencies (Hz)\n");
for iMode = 1:min(maxModesToShow,length(freqHzSorted))
    fprintf("Mode %2d : %12.8e\n", iMode, freqHzSorted(iMode));
end

%% export results to table
outputDir = fileparts(mfilename('fullpath'));
if isempty(outputDir)
    outputDir = pwd;
end
outputDir = fullfile(outputDir,'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% displacement table
dofIndex = (1:GDof)';
nodeIndex = ceil(dofIndex/6);
compIndex = mod(dofIndex-1,6)+1;
compNames = {'UX','UY','UZ','RX','RY','RZ'};
dofComp = compNames(compIndex)';
dispTable = table(dofIndex,nodeIndex,dofComp,displacements, ...
    'VariableNames', {'DOF','Node','Component','Displacement'});
writetable(dispTable, fullfile(outputDir,'ThreeDimFrame_20story_displacements.csv'));

% reactions table
supportDof = prescribedDof;
supportDof = supportDof(:);
reactionNode = ceil(supportDof/6);
reactionComp = compNames(mod(supportDof-1,6)+1)';
reactionSupport = reactionSupport(:);
reactionNode = reactionNode(:);
reactionComp = reactionComp(:);
reactionTable = table(supportDof,reactionNode,reactionComp,reactionSupport, ...
    'VariableNames', {'DOF','Node','Component','Reaction'});
writetable(reactionTable, fullfile(outputDir,'ThreeDimFrame_20story_reactions.csv'));

% modal table
omegaSorted = 2*pi*freqHzSorted;
modeId = (1:numel(freqHzSorted))';
modalTable = table(modeId,freqHzSorted,omegaSorted,sortIdx, ...
    'VariableNames', {'Mode','Frequency_Hz','Omega_rad_s','OriginalIndex'});
writetable(modalTable, fullfile(outputDir,'ThreeDimFrame_20story_modes.csv'));

fprintf('\nResult tables saved to: %s\n', outputDir);
fprintf('- ThreeDimFrame_20story_displacements.csv\n');
fprintf('- ThreeDimFrame_20story_reactions.csv\n');
fprintf('- ThreeDimFrame_20story_modes.csv\n');

%% plotting
plotScale = 40;  % deformation magnification for visibility

% nodal coordinates
nodeDisplacement = [displacements(1:6:end), displacements(2:6:end), ...
    displacements(3:6:end)];
deformedCoords = nodeCoordinates + plotScale * nodeDisplacement;

% undeformed/deformed frame shape
figure('Name','20-story 3D frame shape');
hold on; grid on; axis equal; view(35,20);
for e = 1:numberElements
    n1 = elementNodes(e,1);
    n2 = elementNodes(e,2);
    p1 = nodeCoordinates(n1,:);
    p2 = nodeCoordinates(n2,:);
    p1d = deformedCoords(n1,:);
    p2d = deformedCoords(n2,:);
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k-', 'LineWidth', 1.0);
    plot3([p1d(1) p2d(1)], [p1d(2) p2d(2)], [p1d(3) p2d(3)], ...
        'r--', 'LineWidth', 1.5);
end
xlabel('X'); ylabel('Y'); zlabel('Z');
legend({'Original','Deformed (scaled)'}, 'Location', 'best');
title('3D frame (top load case)');

% story-level drift / top-node displacement by floor
floorLevels = unique(nodeCoordinates(:,2),'stable');
storyUx = zeros(numel(floorLevels),1);
storyUy = zeros(numel(floorLevels),1);
storyUz = zeros(numel(floorLevels),1);
for iFloor = 1:numel(floorLevels)
    floorNodes = find(abs(nodeCoordinates(:,2)-floorLevels(iFloor)) < 1e-9);
    storyUx(iFloor) = mean(abs(nodeDisplacement(floorNodes,1)));
    storyUy(iFloor) = mean(abs(nodeDisplacement(floorNodes,2)));
    storyUz(iFloor) = mean(abs(nodeDisplacement(floorNodes,3)));
end

figure('Name','Story deformation profile');
subplot(1,3,1);
plot(storyUx, floorLevels, 'o-', 'LineWidth', 1.2);
grid on; xlabel('Mean |UX| [m]'); ylabel('Height [m]');
title('Story mean UX');
subplot(1,3,2);
plot(storyUy, floorLevels, 'o-', 'LineWidth', 1.2);
grid on; xlabel('Mean |UY| [m]'); ylabel('Height [m]');
title('Story mean UY');
subplot(1,3,3);
plot(storyUz, floorLevels, 'o-', 'LineWidth', 1.2);
grid on; xlabel('Mean |UZ| [m]'); ylabel('Height [m]');
title('Story mean UZ');

% 1st mode shape in X direction by floor
if ~isempty(modes)
    mode1 = modes(:,1);
    mode1Coords = [mode1(1:6:end), mode1(2:6:end), mode1(3:6:end)];
    firstModeFloor = zeros(numel(floorLevels),1);
    for iFloor = 1:numel(floorLevels)
        floorNodes = find(abs(nodeCoordinates(:,2)-floorLevels(iFloor)) < 1e-9);
        firstModeFloor(iFloor) = mean(abs(mode1Coords(floorNodes,1)));
    end
    if max(firstModeFloor) > 0
        firstModeFloor = firstModeFloor / max(firstModeFloor);
    end
    figure('Name','1st mode shape (normalized, UX component)');
    plot(firstModeFloor, floorLevels, 's-', 'LineWidth', 1.5, 'MarkerSize',6);
    grid on; xlabel('Normalized UX mode shape'); ylabel('Height [m]');
    title('First modal shape by floor');
end

%% local functions
function  [stiffness]=formStiffness3Dframe(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,E,A,Iz,Iy,G,J);
stiffness = zeros(GDof);
for e = 1:numberElements
    indice = elementNodes(e,:);
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

function  [mass] = formMass3Dframe(GDof,numberElements, ...
    elementNodes,numberNodes,nodeCoordinates,rho,A,Iz,Iy)
mass = zeros(GDof);
for e = 1:numberElements
    indice = elementNodes(e,:);
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
    p = (Iz+Iy)/A;
    m = rho*A*L/420*[140 0 0 0 0 0 70 0 0 0 0 0;
                     0  156 0 0 0 22*L 0 54 0 0 0 -13*L;
                     0 0 156 0 -22*L 0 0 0 54 0 13*L 0;
                     0 0 0 140*p 0 0 0 0 0 70*p 0 0;
                     0 0 -22*L 0 4*L^2 0 0 0 -13*L 0 -3*L^2 0;
                     0 22*L 0 0 0 4*L^2 0 13*L 0 0 0 -3*L^2;
                     70 0 0 0 0 0 140 0 0 0 0 0;
                     0 54 0 0 0 13*L 0 156 0 0 0 -22*L;
                     0 0 54 0 -13*L 0 0 0 156 0 22*L 0;
                     0 0 0 70*p 0 0 0 0 0 140*p 0 0;
                     0 0 13*L 0 -3*L^2 0 0 0 22*L 0 4*L^2 0;
                     0 -13*L 0 0 0 -3*L^2 0 -22*L 0 0 0 4*L^2];
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
    mass(elementDof,elementDof)=...
        mass(elementDof,elementDof)+R'*m*R;
end
end

function [Kc,Mc,Fc,T,dofToReduced] = applyRigidDiaphragmConstraints(K, M, F, ...
    nodeCoordinates, startFloor, numFloors)
GDof = size(K,1);
dofToRep = (1:GDof)';
dofTie = [1 3 5]; % UX, UZ, RY (dx, dz, ry)

levels = unique(nodeCoordinates(:,2),'stable');
if numel(levels) < numFloors
    numFloors = numel(levels);
end
if startFloor < 1
    startFloor = 1;
end

for iFloor = startFloor:numFloors
    floorNodes = find(abs(nodeCoordinates(:,2)-levels(iFloor)) < 1e-9);
    if isempty(floorNodes)
        continue;
    end
    repNode = floorNodes(1);
    for dof = dofTie
        repDof = 6*(repNode-1)+dof;
        for j = 1:numel(floorNodes)
            nodeId = floorNodes(j);
            fullDof = 6*(nodeId-1)+dof;
            dofToRep(fullDof) = repDof;
        end
    end
end

uniqueRep = unique(dofToRep,'stable');
dofToReduced = zeros(GDof,1);
for i = 1:numel(uniqueRep)
    dofToReduced(uniqueRep(i)) = i;
end
dofToReduced = dofToReduced(dofToRep);

T = sparse((1:GDof)', dofToReduced, ones(GDof,1), GDof, numel(uniqueRep));
Kc = T' * K * T;
Mc = T' * M * T;
Fc = T' * F;
end

function displacements=solution(GDof,prescribedDof,stiffness,force)
activeDof=setdiff((1:GDof)', [prescribedDof]);
U = stiffness(activeDof,activeDof)\force(activeDof);
displacements = zeros(GDof,1);
displacements(activeDof)=U;
end

function [modes,eigenvalues] = eigenvalue(GDof,prescribedDof, ...
    stiffness,mass,maxEigenvalues)
activeDof = setdiff((1:GDof)', prescribedDof);
if maxEigenvalues == 0
    [V,D] = eig(stiffness(activeDof,activeDof),...
        mass(activeDof,activeDof));
else
    [V,D] = eigs(stiffness(activeDof,activeDof),...
        mass(activeDof,activeDof),maxEigenvalues,'smallestabs');
end
eigenvalues = diag(D);
modes = zeros(GDof,length(eigenvalues));
modes(activeDof,:) = V;
end
