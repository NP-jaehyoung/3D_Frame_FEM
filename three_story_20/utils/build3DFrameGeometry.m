function [nodeCoordinates, elementNodes] = build3DFrameGeometry(numFloors, floorHeight, spanX, spanZ)
% Build regular 20-story 3D frame mesh (square plan with 4 nodes per floor).

if numFloors < 1
    error("numFloors must be at least 1.");
end

nodeCoordinates = zeros(4*numFloors,3);
for i = 1:numFloors
    baseY = floorHeight*(i-1);
    nodeCoordinates(4*(i-1)+1,:) = [0      baseY 0];
    nodeCoordinates(4*(i-1)+2,:) = [0      baseY spanZ];
    nodeCoordinates(4*(i-1)+3,:) = [spanX  baseY spanZ];
    nodeCoordinates(4*(i-1)+4,:) = [spanX  baseY 0];
end

% first floor and middle floors: columns + girders
elementNodes = zeros(8*(numFloors-1),2);
for i = 1:max(numFloors-1,1)
    base = 4*(i-1);
    elementNodes(8*(i-1)+1,:) = [base+1 base+5];
    elementNodes(8*(i-1)+2,:) = [base+2 base+6];
    elementNodes(8*(i-1)+3,:) = [base+3 base+7];
    elementNodes(8*(i-1)+4,:) = [base+4 base+8];
    elementNodes(8*(i-1)+5,:) = [base+5 base+6];
    elementNodes(8*(i-1)+6,:) = [base+6 base+7];
    elementNodes(8*(i-1)+7,:) = [base+7 base+8];
    elementNodes(8*(i-1)+8,:) = [base+8 base+5];
end
end
