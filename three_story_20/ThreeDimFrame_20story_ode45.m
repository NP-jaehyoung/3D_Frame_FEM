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

%% mesh and system matrices
[nodeCoordinates, elementNodes] = build3DFrameGeometry(numFloors, floorHeight, spanX, spanZ);

[stiffness, mass, force, GDof] = ...
    assemble3DFrameMatrices(numFloors, nodeCoordinates, elementNodes, ...
    E, A, Iz, Iy, G, J, rho, 0);
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
