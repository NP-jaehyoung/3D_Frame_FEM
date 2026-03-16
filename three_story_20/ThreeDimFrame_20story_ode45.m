%................................................................
% MATLAB code for 20-story 3D frame dynamic response
% using ode45 with ground acceleration input
clear all
clc
close all

%% add local paths
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, "utils"));

%% model settings
numFloors = 20;            % number of floor levels
floorHeight = 4.0;         % story height [m]
spanX = 4.0;               % x direction span [m]
spanZ = 4.0;               % z direction span [m]
numBaysX = 3;              % number of bays in x
numBaysZ = 3;              % number of bays in z
E = 2e7;                   % modulus of elasticity [kPa]
rho = 7850;                % density
width = 1;                 
depth = 1;                 
nu = 0.167;
[A, Iy, Iz, J, G] = frameSectionProperties(width, depth, E, nu);
floorLoads = zeros(numFloors,1);            % per-floor lateral loads [kN] for static preload
floorLoadDistribution = "AREA";             % UNIFORM or AREA
% floorLoads = -3000 * ones(numFloors,1);   % (uncomment for distributed static lateral loads)
includeDeadLoadAsMass = true;              % convert dead load to mass contribution
floorDeadLoadKiloNewton = zeros(numFloors,1); % dead load [kN] per floor, positive for downward load
deadLoadDistribution = "AREA";             % UNIFORM or AREA
deadLoadAcceleration = 9.81;               % conversion g (m/s^2)
includeStaticForcesInDynamic = false;       % keep ground-motion-only response in dynamic script

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

%% mesh and system matrices
[nodeCoordinates, elementNodes, floorNodeIds] = build3DFrameGeometry( ...
    numFloors, floorHeight, spanX, spanZ, numBaysX, numBaysZ);

[stiffness, mass, force, GDof, topNode] = ...
    assemble3DFrameMatrices(numFloors, nodeCoordinates, elementNodes, ...
    E, A, Iz, Iy, G, J, rho, 0, floorLoads, "UX", floorNodeIds, floorLoadDistribution);

if includeDeadLoadAsMass
    floorMass = abs(floorDeadLoadKiloNewton) / deadLoadAcceleration;
    massAll = addDeadLoadMassToLumpedNodes(numFloors, mass, nodeCoordinates, ...
        floorNodeIds, floorMass, deadLoadDistribution);
else
    massAll = mass;
end
stiffnessAll = stiffness;
forceAll = force;

if includeDeadLoadAsMass
    totalDeadMass = sum(floorMass);
    fprintf("Dead load converted to added mass: total = %.6f [mass unit], distribution = %s\n", ...
        totalDeadMass, deadLoadDistribution);
end

%% rigid diaphragm constraints: tie each floor as a rigid body about its centroid
[stiffness, mass, force, T, dofToReduced] = ...
    applyRigidDiaphragmConstraints(stiffnessAll, massAll, forceAll, ...
    nodeCoordinates, 2, numFloors);

%% fixed support: all base-floor nodes
baseNodes = floorNodeIds{1}(:);
prescribedDof = reshape((6*(baseNodes - 1) + (1:6))', 1, []);
prescribedDofReduced = unique(dofToReduced(prescribedDof));
activeDof = setdiff((1:size(stiffness,1))', prescribedDofReduced);

% reduced active matrices
Kff = stiffness(activeDof, activeDof);
Mff = mass(activeDof, activeDof);

%% modal-based Rayleigh damping (using reduced matrices)
[~, eigenvalues] = solveEigenModes(size(stiffness,1), prescribedDofReduced, ...
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
odeRhs = @(t, y) ode45StateRhs3DFrame(t, y, Mff, Cff, Kff, baseLoadActive, agFun);
y0 = zeros(2*size(Kff,1),1);    % [u; v]

if includeStaticForcesInDynamic && any(abs(forceAll) > 1e-12)
    forceReduced = T' * forceAll;
    forceActive = forceReduced(activeDof);
    y0(1:size(Kff,1),1) = Mff \ forceActive;
end
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
if includeStaticForcesInDynamic
    fprintf("Dynamic response includes static force preload (`floorLoads`, `topLoad`)\n");
else
    fprintf("Dynamic response uses ground acceleration only (no static preload force)\n");
end

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

function massAll = addDeadLoadMassToLumpedNodes(numFloors, massAll, nodeCoordinates, ...
    floorNodeIds, floorMass, distribution)
if nargin < 6 || isempty(distribution)
    distribution = "UNIFORM";
end
distribution = upper(string(distribution));
if isscalar(floorMass)
    floorMass = floorMass * ones(numFloors,1);
end
if numel(floorMass) ~= numFloors
    error("addDeadLoadMassToLumpedNodes: floorMass size must match numFloors");
end
if ischar(floorMass) || isstring(floorMass)
    error("addDeadLoadMassToLumpedNodes: floorMass must be numeric.");
end

switch distribution
    case {"UNIFORM","EQUAL","POINT"}
        distribution = "UNIFORM";
    case {"AREA","GEOMETRIC","GEOMETRY"}
        distribution = "AREA";
    otherwise
        error("addDeadLoadMassToLumpedNodes: unsupported distribution '%s'", distribution);
end

if isempty(floorNodeIds)
    error("addDeadLoadMassToLumpedNodes: floorNodeIds is empty.");
end

for iFloor = 1:numFloors
    if isempty(floorNodeIds{iFloor})
        continue;
    end
    floorNodes = floorNodeIds{iFloor};
    if distribution == "UNIFORM"
        nodeWeights = ones(numel(floorNodes),1) / numel(floorNodes);
    else
        nodeWeights = computeAreaLikeWeights(nodeCoordinates(floorNodes, [1,3]));
    end
    if abs(sum(nodeWeights) - 1) > 1e-10
        nodeWeights = nodeWeights / max(sum(nodeWeights),eps);
    end
    nodeMass = floorMass(iFloor) * nodeWeights(:);
    for j = 1:numel(floorNodes)
        dofs = 6*(floorNodes(j)-1) + (1:3);
        massAll(dofs, dofs) = massAll(dofs, dofs) + diag(nodeMass(j) * ones(3,1));
    end
end
end

function nodeWeights = computeAreaLikeWeights(nodeXZ)
numNodes = size(nodeXZ,1);
if numNodes == 0
    nodeWeights = zeros(0,1);
    return;
end
if numNodes <= 2
    nodeWeights = ones(numNodes,1) / numNodes;
    return;
end

xyAll = nodeXZ(:,1:2);
nodeIds = (1:numNodes)';

theta0 = atan2(xyAll(:,2) - mean(xyAll(:,2)), xyAll(:,1) - mean(xyAll(:,1)));
[~, order0] = sort(theta0);
xy = xyAll(order0,:);
nodeIds = nodeIds(order0);

area2 = 0.0;
for i = 1:numNodes
    j = mod(i,numNodes) + 1;
    area2 = area2 + xy(i,1)*xy(j,2) - xy(j,1)*xy(i,2);
end
if abs(area2) < 1e-14
    nodeWeights = ones(numNodes,1) / numNodes;
    return;
end

if area2 < 0
    xy = flipud(xy);
    area2 = -area2;
    nodeIds = flipud(nodeIds);
end

centroid = [0, 0];
for i = 1:numNodes
    j = mod(i,numNodes) + 1;
    crossVal = xy(i,1)*xy(j,2) - xy(j,1)*xy(i,2);
    centroid(1) = centroid(1) + (xy(i,1) + xy(j,1)) * crossVal;
    centroid(2) = centroid(2) + (xy(i,2) + xy(j,2)) * crossVal;
end
centroid = centroid / (3*area2);

theta = atan2(xy(:,2)-centroid(2), xy(:,1)-centroid(1));
[~, order] = sort(theta);
xy = xy(order,:);
nodeIds = nodeIds(order);

edgeSector = zeros(numNodes,1);
for i = 1:numNodes
    j = mod(i,numNodes) + 1;
    p = xy(i,:) - centroid;
    q = xy(j,:) - centroid;
    triArea = 0.5 * abs(p(1)*q(2) - p(2)*q(1));
    edgeSector(i) = edgeSector(i) + 0.5*triArea;
    edgeSector(j) = edgeSector(j) + 0.5*triArea;
end

if all(edgeSector < 1e-14)
    nodeWeights = ones(numNodes,1) / numNodes;
    return;
end

nodeWeightsSorted = edgeSector / sum(edgeSector);

% Restore original node ordering
restoreNodeIds = nodeIds;
nodeWeights = zeros(numNodes,1);
nodeWeights(restoreNodeIds) = nodeWeightsSorted;
end
