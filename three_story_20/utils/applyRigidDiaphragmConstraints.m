function [stiffness, mass, force, T, dofToReduced] = ...
    applyRigidDiaphragmConstraints(stiffness, mass, force, ...
    nodeCoordinates, startFloor, numFloors)

GDof = size(stiffness,1);
dofToRep = (1:GDof)';
dofTie = [1 3 5]; % UX, UZ, RY (dx, dz, ry)

levels = unique(nodeCoordinates(:,2),'stable');
if numel(levels) < numFloors
    numFloors = numel(levels);
end
if nargin < 5 || isempty(startFloor)
    startFloor = 1;
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
stiffness = T' * stiffness * T;
mass = T' * mass * T;
force = T' * force;
end
