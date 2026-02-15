%................................................................
% MATLAB code for 20-story 3D frame dynamic response
% using ode45 with ground acceleration input
clear all
clc
close all

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

%% time history settings
timeEnd = 20;               % sec
dt = 0.01;                  % integration step for target output
time = (0:dt:timeEnd)';
% sample ground acceleration in m/s^2
agPeak = 0.30*9.81;         % 0.30g
agFreq = 1.8;               % Hz
agPulseCenter = 5.0;        % s
agPulseWidth = 1.1;         % s
groundAcceleration = agPeak .* exp(-((time-agPulseCenter)/agPulseWidth).^2) ...
    .* sin(2*pi*agFreq*time);

%% damping settings
dampingRatio = 0.02;        % target damping ratio (ζ)
zeta1 = dampingRatio;

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
force = zeros(GDof,1);  % no external active load in this example
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

%% fixed support: nodes 1~4 on first floor
prescribedDof = 1:24;
prescribedDofReduced = unique(dofToReduced(prescribedDof));
activeDof = setdiff((1:size(stiffness,1))', prescribedDofReduced);

% reduced active matrices
Kff = stiffness(activeDof, activeDof);
Mff = mass(activeDof, activeDof);

%% modal-based Rayleigh damping (using reduced matrices)
[~, eigenvalues] = eigenvalue(size(stiffness,1), prescribedDofReduced, ...
    stiffness, mass, 6);
omegaAll = sqrt(abs(eigenvalues));
omegaAll = sort(omegaAll(omegaAll > 1e-8));
if numel(omegaAll) >= 2
    w1 = omegaAll(1);
    w2 = omegaAll(2);
    alphaM = 2*zeta1*w1*w2/(w1 + w2);
    betaK = 2*zeta1/(w1 + w2);
elseif numel(omegaAll) == 1
    w1 = omegaAll(1);
    alphaM = 0;
    betaK = 2*zeta1/w1;
else
    error("No natural frequency found for damping calibration.");
end
Cfull = alphaM*mass + betaK*stiffness;
Cff = Cfull(activeDof, activeDof);

%% ground excitation vector (all x-translation DOFs = 1)
groundDirectionFull = zeros(GDof,1);
groundDirectionFull(1:6:end) = 1;   % UX of every node
baseLoadReduced = T' * (-massAll * groundDirectionFull);
baseLoadActive = baseLoadReduced(activeDof);

%% ode45 integration
agFun = @(t) interp1(time, groundAcceleration, t, 'pchip', 0);
odeRhs = @(t, y) dynamicStateRhs(t, y, Mff, Cff, Kff, baseLoadActive, agFun);
y0 = zeros(2*size(Kff,1),1);    % [u; v]
odeOptions = odeset('RelTol',1e-6, 'AbsTol',1e-8);
[timeSol, ySol] = ode45(odeRhs, time, y0, odeOptions);

uAct = ySol(:,1:size(Kff,1));         % nActive x nTime (rows)
vAct = ySol(:,size(Kff,1)+1:end);     % nActive x nTime (rows)
agSol = agFun(timeSol);

% reduced-space acceleration
nSteps = numel(timeSol);
accActive = zeros(nSteps, size(Kff,1));
for iTime = 1:nSteps
    rhsForce = baseLoadActive * agSol(iTime);
    accActive(iTime,:) = (Mff \ (-Cff*vAct(iTime,:)' - Kff*uAct(iTime,:)' + rhsForce))';
end

%% reconstruct full (reduced-space and full-space)
nReduced = size(stiffness,1);
uReduced = zeros(nSteps, nReduced);
vReduced = zeros(nSteps, nReduced);
aReduced = zeros(nSteps, nReduced);
uReduced(:,activeDof) = uAct;
vReduced(:,activeDof) = vAct;
aReduced(:,activeDof) = accActive;

displacements = (T * uReduced')';
velocities    = (T * vReduced')';
accelerations = (T * aReduced')';

%% response extraction
topNode = 4*(numFloors-1) + 3;   % 1st corner node on top floor
topUx = displacements(:,6*topNode-5);
topUy = displacements(:,6*topNode-4);
topUz = displacements(:,6*topNode-3);
topAx = accelerations(:,6*topNode-5);
topAy = accelerations(:,6*topNode-4);
topAz = accelerations(:,6*topNode-3);

levels = unique(nodeCoordinates(:,2),'stable');
nLevels = numel(levels);
storyDrift = zeros(nLevels-1, nSteps);
for iFloor = 2:nLevels
    floorNodesTop = find(abs(nodeCoordinates(:,2)-levels(iFloor)) < 1e-9);
    floorNodesBottom = find(abs(nodeCoordinates(:,2)-levels(iFloor-1)) < 1e-9);
    uxTop = displacements(:, 6*(floorNodesTop)-5);
    uxBottom = displacements(:, 6*(floorNodesBottom)-5);
    storyDrift(iFloor-1,:) = mean(uxTop - uxBottom,2);
end
maxStoryDrift = max(abs(storyDrift),[],2);
[maxDriftValue, maxDriftFloorIdx] = max(maxStoryDrift);
maxDriftFloor = levels(maxDriftFloorIdx+1);

fprintf("\n20-story dynamic response (ode45) example\n");
fprintf("Solved time points: %d (%.3f ~ %.3f) sec\n", ...
    nSteps, timeSol(1), timeSol(end));
fprintf("Peak top UX = %.6e m\n", max(abs(topUx)));
fprintf("Peak top UY = %.6e m\n", max(abs(topUy)));
fprintf("Peak top UZ = %.6e m\n", max(abs(topUz)));
fprintf("Max story drift=%.6e m at story %.1f m\n", maxDriftValue, maxDriftFloor);

%% plotting
figure('Name','Input ground acceleration');
plot(timeSol, agSol/9.81, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Ground Acceleration [g]');
title('Input ground acceleration');

figure('Name','Top node displacement history');
plot(timeSol, topUx, 'LineWidth', 1.2); hold on;
plot(timeSol, topUy, 'LineWidth', 1.0);
plot(timeSol, topUz, 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Displacement [m]');
legend({'UX','UY','UZ'}, 'Location', 'best');
title(sprintf('Top node (node %d) displacements', topNode));

figure('Name','Top node acceleration history');
plot(timeSol, topAx/9.81, 'LineWidth', 1.2); hold on;
plot(timeSol, topAy/9.81, 'LineWidth', 1.0);
plot(timeSol, topAz/9.81, 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Acceleration [g]');
legend({'AX','AY','AZ'}, 'Location', 'best');
title(sprintf('Top node (node %d) accelerations', topNode));

figure('Name','Story drift histories');
plot(timeSol, storyDrift', 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Inter-story drift (UX) [m]');
legend(arrayfun(@(k) sprintf('Story %d', k), 2:nLevels, 'UniformOutput', false), ...
    'Location', 'eastoutside');
title('Story drift in UX');

%% export results to table
outputDir = fileparts(mfilename('fullpath'));
if isempty(outputDir)
    outputDir = pwd;
end
outputDir = fullfile(outputDir,'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

topTable = table(timeSol, agSol, topUx, topUy, topUz, topAx, topAy, topAz, ...
    'VariableNames', {'Time_s','GroundAcc_m_per_s2','TopUx_m','TopUy_m','TopUz_m',...
    'TopAx_m_per_s2','TopAy_m_per_s2','TopAz_m_per_s2'});
writetable(topTable, fullfile(outputDir,'ThreeDimFrame_20story_ode45_topHistory.csv'));

storyNames = arrayfun(@(k) sprintf('Story%02d',k), 2:nLevels, 'UniformOutput', false);
storyTable = array2table([timeSol, storyDrift'], ...
    'VariableNames', [{'Time_s'}, storyNames]);
writetable(storyTable, fullfile(outputDir,'ThreeDimFrame_20story_ode45_storyDrift.csv'));

fprintf('\nResult tables saved to: %s\n', outputDir);
fprintf('- ThreeDimFrame_20story_ode45_topHistory.csv\n');
fprintf('- ThreeDimFrame_20story_ode45_storyDrift.csv\n');

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

function dYdt = dynamicStateRhs(t, y, M, C, K, groundLoad, agFun)
n = size(M,1);
u = y(1:n);
v = y(n+1:end);
ag = agFun(t);
rhs = -C*v - K*u + groundLoad*ag;
a = M \ rhs;
dYdt = [v; a];
end
