function [nodeCoordinates, elementNodes, floorNodeIds, topCornerNode] = build3DFrameGeometry(numFloors, floorHeight, spanX, spanZ, numBaysX, numBaysZ)
% Build a regular multi-bay 3D frame mesh on a rectangular plan grid.

if numFloors < 2
    error("numFloors must be at least 2 to form a frame.");
end

if nargin < 5 || isempty(numBaysX)
    numBaysX = 3;
end
if nargin < 6 || isempty(numBaysZ)
    numBaysZ = 3;
end
if numBaysX < 1 || numBaysZ < 1
    error("numBaysX and numBaysZ must be at least 1.");
end

nodesPerLevel = (numBaysX + 1) * (numBaysZ + 1);
nodeCoordinates = zeros(nodesPerLevel * numFloors, 3);
floorNodeIds = cell(numFloors, 1);

nodeId = 0;
for iFloor = 1:numFloors
    baseY = floorHeight * (iFloor - 1);
    localIds = zeros(1, nodesPerLevel);
    localIdx = 0;
    for iz = 0:numBaysZ
        for ix = 0:numBaysX
            nodeId = nodeId + 1;
            localIdx = localIdx + 1;
            nodeCoordinates(nodeId, :) = [ix * spanX, baseY, iz * spanZ];
            localIds(localIdx) = nodeId;
        end
    end
    floorNodeIds{iFloor} = localIds;
end

% Roof corner at max X and Z = 0, matching the previous single-bay loading corner.
topCornerNode = floorNodeIds{end}(numBaysX + 1);

elementCount = (numFloors - 1) * ...
    (nodesPerLevel + (numBaysZ + 1) * numBaysX + (numBaysX + 1) * numBaysZ);
elementNodes = zeros(elementCount, 2);
eIdx = 0;

for iFloor = 1:(numFloors - 1)
    floorBottom = floorNodeIds{iFloor};
    floorTop = floorNodeIds{iFloor + 1};

    for localIdx = 1:nodesPerLevel
        eIdx = eIdx + 1;
        elementNodes(eIdx, :) = [floorBottom(localIdx), floorTop(localIdx)];
    end

    for iz = 0:numBaysZ
        row = iz * (numBaysX + 1);
        for ix = 1:numBaysX
            eIdx = eIdx + 1;
            n1 = floorTop(row + ix);
            n2 = floorTop(row + ix + 1);
            elementNodes(eIdx, :) = [n1, n2];
        end
    end

    for iz = 0:(numBaysZ - 1)
        row = iz * (numBaysX + 1);
        nextRow = (iz + 1) * (numBaysX + 1);
        for ix = 1:(numBaysX + 1)
            eIdx = eIdx + 1;
            n1 = floorTop(row + ix);
            n2 = floorTop(nextRow + ix);
            elementNodes(eIdx, :) = [n1, n2];
        end
    end
end

elementNodes = elementNodes(1:eIdx, :);
end
