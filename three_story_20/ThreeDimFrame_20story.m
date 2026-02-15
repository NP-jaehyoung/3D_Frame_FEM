%................................................................
% MATLAB codes for Finite Element Analysis
% 20-story 3D frame example
% clear memory
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
if width <= 0 || depth <= 0
    error("width and depth must be positive.");
end
Iy = width*depth^3/12;      % second moment about local y
Iz = depth*width^3/12;      % second moment about local z
nu = 0.167;
if width >= depth
    ar = depth/width;
    J = depth * width^3 * (1/3 - 0.21*ar + 0.063*ar^5);
else
    ar = width/depth;
    J = width * depth^3 * (1/3 - 0.21*ar + 0.063*ar^5);
end
if J <= 0
    error("Invalid torsional constant J (non-positive). Check width/depth.");
end
G = E/(2*(1+nu));
floorLoads = -3000 * ones(numFloors,1);   % floor lateral load per level [kN]
topLoad = floorLoads(end);               % for log only: top-floor equivalent point load

%% mesh and loads
[nodeCoordinates, elementNodes] = build3DFrameGeometry(numFloors, floorHeight, spanX, spanZ);
numberElements = size(elementNodes,1);
topNode = 4*(numFloors-1) + 3;
topLoadDof = 6*topNode - 5;

%% assemble K, M
[stiffness, mass, force, GDof] = ...
    assemble3DFrameMatrices(numFloors, nodeCoordinates, elementNodes, ...
    E, A, Iz, Iy, G, J, rho, 0, floorLoads, "UX");
stiffnessAll = stiffness;
massAll = mass;
forceAll = force;

%% rigid diaphragm constraints: full rigid-floor coupling (UX,UY,UZ with RX,RY,RZ)
[stiffness, mass, force, T, dofToReduced] = ...
    applyRigidDiaphragmConstraints(stiffnessAll, massAll, forceAll, ...
    nodeCoordinates, 2, numFloors);

%% boundary conditions: fix 1st floor nodes (nodes 1~4)
prescribedDof = 1:24;
prescribedDofReduced = unique(dofToReduced(prescribedDof));

%% static solution and output
displacementsReduced = solveConstrainedStatic(size(stiffness,1), ...
    prescribedDofReduced, ...
    stiffness, force);
displacements = T * displacementsReduced;
reactions = stiffnessAll * displacements - forceAll;
reactionSupport = reactions(prescribedDof);

fprintf("20-floor 3D frame example\n");
fprintf("Top node: %d, top load DOF: %d, load: %.4g\n", topNode, topLoadDof, topLoad);
fprintf("Floor loads: applied per floor (each floor equally distributed to 4 nodes)\n");
fprintf("Top node displacement (UX,UY,UZ) [m]\n");
fprintf("UX = %.6e\n", displacements(topLoadDof));
fprintf("UY = %.6e\n", displacements(topLoadDof+1));
fprintf("UZ = %.6e\n\n", displacements(topLoadDof+2));

fprintf("Support reactions (first 24 DOF)\n");
for ii = 1:numel(prescribedDof)
    fprintf("DOF %d : %+13.6e\n", prescribedDof(ii), reactionSupport(ii));
end

%% support reaction summary: translational load extraction prepared
fprintf("Applied load extraction will be written after GDof definition.\n");

%% modal analysis
maxModesToShow = 10;
[modesReduced,eigenvalues] = solveEigenModes(size(stiffness,1), ...
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

% floor load table (user input summary)
floorTable = table((1:numFloors)', floorLoads, ...
    'VariableNames', {'Floor','FloorLoad_kN'});
writetable(floorTable, fullfile(outputDir,'ThreeDimFrame_20story_floorLoads.csv'));

% nodal rotation table (rigid floor rotations)
ryNodeIds = (1:size(nodeCoordinates,1))';
ryDof = 6*ryNodeIds - 1;
rzDof = 6*ryNodeIds;
rxDof = 6*ryNodeIds - 2;
rxValue = displacements(rxDof);
ryValue = displacements(ryDof);
rzValue = displacements(rzDof);
rotTable = table(rxDof, ryDof, rzDof, ryNodeIds, rxValue, ryValue, rzValue, ...
    'VariableNames', {'DOF_RX','DOF_RY','DOF_RZ','Node','RX_Rad','RY_Rad','RZ_Rad'});
writetable(rotTable, fullfile(outputDir,'ThreeDimFrame_20story_node_rotations.csv'));

% applied load table (per DOF)
loadDof = find(abs(forceAll) > 1e-12);
loadNode = ceil(loadDof/6);
loadCompIndex = mod(loadDof-1,6)+1;
loadValue = forceAll(loadDof);
loadCompName = compNames(loadCompIndex)';

loadDof = loadDof(:);
loadNode = loadNode(:);
loadCompIndex = loadCompIndex(:);
loadValue = loadValue(:);
loadCompName = loadCompName(:);
loadDirection = loadCompName(loadCompIndex <= 3); % translational for arrows
loadMoment = loadCompName(loadCompIndex >= 4);    % rotational DOFs (shown in table only)

loadTable = table(loadDof,loadNode,loadCompName,loadValue, ...
    'VariableNames', {'DOF','Node','Component','AppliedLoad'});
writetable(loadTable, fullfile(outputDir,'ThreeDimFrame_20story_external_loads.csv'));

fprintf("Applied loads (non-zero DOF)\n");
for ii = 1:numel(loadDof)
    fprintf("DOF %d (Node %d, %s) : %+13.6e [kN]\n", ...
        loadDof(ii), loadNode(ii), loadCompName{ii}, loadValue(ii));
end

fprintf('\nResult tables saved to: %s\n', outputDir);
fprintf('- ThreeDimFrame_20story_displacements.csv\n');
fprintf('- ThreeDimFrame_20story_reactions.csv\n');
fprintf('- ThreeDimFrame_20story_modes.csv\n');
fprintf('- ThreeDimFrame_20story_floorLoads.csv\n');
fprintf('- ThreeDimFrame_20story_external_loads.csv\n');
fprintf('- ThreeDimFrame_20story_node_rotations.csv\n');

%% plotting
plotScale = 40;  % deformation magnification for visibility

floorLevels = unique(nodeCoordinates(:,2),'stable');

% nodal coordinates
nodeDisplacement = [displacements(1:6:end), displacements(2:6:end), ...
    displacements(3:6:end)];
deformedCoords = nodeCoordinates + plotScale * nodeDisplacement;
% torsion-aware coordinates (add rigid floor RY rotation effect)
floorTorsionCoords = zeros(size(nodeCoordinates));
for iFloor = 1:numel(floorLevels)
    floorNodes = find(abs(nodeCoordinates(:,2)-floorLevels(iFloor)) < 1e-9);
    if isempty(floorNodes)
        continue;
    end
    floorCentroid = mean(nodeCoordinates(floorNodes,:),1);
    theta = mean(displacements(6*floorNodes-1)); % representative floor RY
    c = cos(theta);
    s = sin(theta);
    for j = 1:numel(floorNodes)
        nId = floorNodes(j);
        basePoint = nodeCoordinates(nId,:);
        localVec = basePoint - floorCentroid;
        transVec = [displacements(6*nId-5), displacements(6*nId-4), displacements(6*nId-3)];
        rotVec = [c*localVec(1) + s*localVec(3), 0, -s*localVec(1) + c*localVec(3)];
        floorTorsionCoords(nId,:) = (floorCentroid + transVec) + [rotVec(1)-localVec(1), transVec(2), rotVec(3)-localVec(3)];
    end
end

% floor polygon connectivity
quadOrder = [1,2,3,4,1];

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

% draw applied nodal loads (translational components only)
coordSpan = max(nodeCoordinates) - min(nodeCoordinates);
baseArrowScale = max(coordSpan);
if isempty(baseArrowScale) || baseArrowScale == 0
    baseArrowScale = 1.0;
end
transLoadIdx = find(loadCompIndex <= 3);
if ~isempty(transLoadIdx)
    maxTransLoad = max(abs(loadValue(transLoadIdx)));
    if maxTransLoad > 0
        arrowScale = 0.30 * baseArrowScale / maxTransLoad;
    else
        arrowScale = 0;
    end
    for ii = transLoadIdx'
        dirVec = zeros(1,3);
        switch loadCompIndex(ii)
            case 1, dirVec = [1,0,0];
            case 2, dirVec = [0,1,0];
            case 3, dirVec = [0,0,1];
        end
        p = nodeCoordinates(loadNode(ii),:);
        q = dirVec * (loadValue(ii) * arrowScale);
        plotTag = "";
        if abs(loadValue(ii)) >= 1e-12
            plotTag = sprintf('%s=%+.1f', loadCompName{ii}, loadValue(ii));
            quiver3(p(1), p(2), p(3), q(1), q(2), q(3), 0, ...
                'Color',[0.05 0.25 0.85], 'LineWidth',1.5, 'MaxHeadSize',1.0);
            text(p(1)+0.5*q(1), p(2)+0.5*q(2), p(3)+0.5*q(3), plotTag, ...
                'FontSize',8, 'Color',[0.05 0.25 0.85], 'FontWeight','bold');
        end
    end
end

xlabel('X'); ylabel('Y'); zlabel('Z');
if isempty(loadDirection)
    legend({'Original','Deformed (scaled)'}, 'Location', 'best');
else
    legend({'Original','Deformed (scaled)', 'Load vector'}, 'Location', 'best');
end
title('3D frame (top load case)');

% torsion shape at each floor (in-plane rigid rotation RY)
figure('Name','20-story floor torsion');
hold on; grid on; axis equal; view(35,25);
for iFloor = 1:numel(floorLevels)
    floorNodes = (4*(iFloor-1)+1):(4*iFloor);
    baseFloor = nodeCoordinates(floorNodes,:);
    torsFloor = floorTorsionCoords(floorNodes,:);
    plot3(baseFloor(quadOrder,1), baseFloor(quadOrder,2), baseFloor(quadOrder,3), ...
        'k-', 'LineWidth',1.0);
    plot3(torsFloor(quadOrder,1), torsFloor(quadOrder,2), torsFloor(quadOrder,3), ...
        'm--', 'LineWidth',2.0);
    c = mean(torsFloor(:,1));
    d = mean(torsFloor(:,2));
    e = mean(torsFloor(:,3));
    text(c, d, e, sprintf('f%02d RY=%+.4e', iFloor, mean(displacements(floorNodes*6-1))), ...
        'FontSize', 8, 'Color',[0.5 0 0.5], 'FontWeight','bold');
end
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Rigid-floor twist (RY) visualization');
legend({'Original floor','RY-rotated floor'}, 'Location', 'best');

% story-level drift / top-node displacement by floor
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
