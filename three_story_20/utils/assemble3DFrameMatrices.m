function [stiffness, mass, force, GDof, topNode, topLoadDof] = ...
    assemble3DFrameMatrices(numFloors, nodeCoordinates, elementNodes, ...
    E, A, Iz, Iy, G, J, rho, topLoad, floorLoad, floorLoadDof, ...
    floorNodeIds, floorLoadDistribution)
numberNodes = size(nodeCoordinates,1);
numberElements = size(elementNodes,1);
GDof = 6*numberNodes;

if nargin < 10 || isempty(topLoad)
    topLoad = 0;
end

if nargin < 11 || isempty(floorLoad)
    floorLoad = zeros(numFloors,1);
end
if isscalar(floorLoad)
    floorLoad = floorLoad * ones(numFloors,1);
end
if numel(floorLoad) ~= numFloors
    error("assemble3DFrameMatrices: floorLoad size must match numFloors");
end

if nargin < 12 || isempty(floorLoadDof)
    floorLoadDof = 1;  % UX
end
if ischar(floorLoadDof) || isstring(floorLoadDof)
    switch upper(string(floorLoadDof))
        case "UX", floorLoadDof = 1;
        case "UY", floorLoadDof = 2;
        case "UZ", floorLoadDof = 3;
        case "RX", floorLoadDof = 4;
        case "RY", floorLoadDof = 5;
        case "RZ", floorLoadDof = 6;
        otherwise, error("assemble3DFrameMatrices: unknown floorLoadDof '%s'", floorLoadDof);
    end
end
if ~isscalar(floorLoadDof) || floorLoadDof < 1 || floorLoadDof > 6
    error("assemble3DFrameMatrices: floorLoadDof must be 1..6");
end

if nargin < 13 || isempty(floorNodeIds)
    levels = unique(nodeCoordinates(:,2),'stable');
    floorNodeIds = cell(numel(levels),1);
    for iFloor = 1:numel(levels)
        floorNodes = find(abs(nodeCoordinates(:,2)-levels(iFloor)) < 1e-9);
        if isempty(floorNodes)
            floorNodeIds{iFloor} = [];
        else
            floorNodeIds{iFloor} = floorNodes(:);
        end
    end
end
if iscell(floorNodeIds)
    if numel(floorNodeIds) < numFloors
        error("assemble3DFrameMatrices: floorNodeIds has fewer floors than numFloors");
    end
    floorNodeIds = floorNodeIds(1:numFloors);
else
    if size(floorNodeIds,1) == 1 && numFloors > 1
        error("assemble3DFrameMatrices: floorNodeIds must be provided per floor");
    end
    temp = cell(numFloors,1);
    if size(floorNodeIds,1) == numFloors
        for iFloor = 1:numFloors
            ids = floorNodeIds(iFloor,:);
            temp{iFloor} = ids(ids ~= 0);
        end
    elseif size(floorNodeIds,1) > numFloors
    for iFloor = 1:numFloors
        temp{iFloor} = floorNodeIds(iFloor,:);
        temp{iFloor} = temp{iFloor}(temp{iFloor} ~= 0);
        end
    else
        error("assemble3DFrameMatrices: floorNodeIds size is not recognized");
    end
    floorNodeIds = temp;
end

if nargin < 14 || isempty(floorLoadDistribution)
    floorLoadDistribution = "UNIFORM";
end
if ischar(floorLoadDistribution) || isstring(floorLoadDistribution)
    floorLoadDistribution = upper(string(floorLoadDistribution));
else
    error("assemble3DFrameMatrices: floorLoadDistribution must be string");
end
switch floorLoadDistribution
    case {"UNIFORM","EQUAL","POINT"}
        floorLoadDistribution = "UNIFORM";
    case {"AREA","GEOMETRIC","GEOMETRY"}
        floorLoadDistribution = "AREA";
    otherwise
        error("assemble3DFrameMatrices: unsupported floorLoadDistribution '%s'", floorLoadDistribution);
end

if isempty(floorNodeIds)
    error("assemble3DFrameMatrices: floorNodeIds is empty.");
end
topFloorNodes = floorNodeIds{min(numFloors, numel(floorNodeIds))};
if isempty(topFloorNodes)
    error("assemble3DFrameMatrices: top floor has no node list.");
end

% top node load (positive in +X; negative value means -X direction)
topNode = topFloorNodes(1);
topLoadDof = 6*topNode - 5;      % UX of topNode
force = zeros(GDof,1);
if topLoad ~= 0
    force(topLoadDof) = topLoad;
end

% floor loads: applied at each floor, distributed to each floor nodes
for iFloor = 1:numFloors
    loadValue = floorLoad(iFloor);
    if loadValue == 0
        continue;
    end

    floorNodes = floorNodeIds{iFloor}(:);
    if isempty(floorNodes)
        error("assemble3DFrameMatrices: no nodes defined for floor %d", iFloor);
    end

    if floorLoadDistribution == "UNIFORM"
        nodeWeights = ones(numel(floorNodes),1) / numel(floorNodes);
    else
        nodeWeights = computeAreaLikeWeights(nodeCoordinates(floorNodes, [1,3]));
    end

    floorDof = 6*(floorNodes - 1) + floorLoadDof;
    force(floorDof) = force(floorDof) + (loadValue * nodeWeights);
end

stiffness = formStiffness3DframeInternal(GDof, numberElements, elementNodes, ...
    nodeCoordinates, E, A, Iz, Iy, G, J);
mass = formMass3DframeInternal(GDof, numberElements, elementNodes, ...
    nodeCoordinates, rho, A, Iz, Iy);
end

function  [stiffness]=formStiffness3DframeInternal(GDof,numberElements,elementNodes,...
    nodeCoordinates,E,A,Iz,Iy,G,J)
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

function  [mass] = formMass3DframeInternal(GDof,numberElements, ...
    elementNodes,nodeCoordinates,rho,A,Iz,Iy)
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

% Sort nodes to form floor polygon in x-z plane and compute centroid
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

% Reorder once more around centroid so sector partition is consistent
theta = atan2(xy(:,2)-centroid(2), xy(:,1)-centroid(1));
[~, order] = sort(theta);
xy = xy(order,:);

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
restoreNodeIds = nodeIds(order);
nodeWeights = zeros(numNodes,1);
nodeWeights(restoreNodeIds) = nodeWeightsSorted;
end
