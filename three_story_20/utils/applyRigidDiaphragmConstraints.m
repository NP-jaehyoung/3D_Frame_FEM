function [stiffness, mass, force, T, dofToReduced] = ...
    applyRigidDiaphragmConstraints(stiffness, mass, force, ...
    nodeCoordinates, startFloor, numFloors)

GDof = size(stiffness,1);

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
if startFloor > numFloors
    startFloor = numFloors;
end

% build reduced-coordinate map with rigid diaphragm equations
% UX_i = UXc - (z_i-zc)*RY
% UZ_i = UZc + (x_i-xc)*RY
% RY_i = RY
rows = [];
cols = [];
vals = [];
colIdx = 0;

for iFloor = startFloor:numFloors
    floorNodes = find(abs(nodeCoordinates(:,2)-levels(iFloor)) < 1e-9);
    if isempty(floorNodes)
        continue;
    end

    % rigid-floor generalized DOFs (centroid translations + twist)
    cx = mean(nodeCoordinates(floorNodes,1));
    cz = mean(nodeCoordinates(floorNodes,3));
    colUx = colIdx + 1; colIdx = colIdx + 1;  % centroid UX
    colUz = colIdx + 1; colIdx = colIdx + 1;  % centroid UZ
    colRy = colIdx + 1; colIdx = colIdx + 1;  % floor twist

    for j = 1:numel(floorNodes)
        nodeId = floorNodes(j);
        xj = nodeCoordinates(nodeId,1);
        zj = nodeCoordinates(nodeId,3);
        dx = xj - cx;
        dz = zj - cz;

        dofUx = 6*(nodeId-1) + 1;
        dofUy = dofUx + 1;
        dofUz = dofUx + 2;
        dofRx = dofUx + 3;
        dofRy = dofUx + 4;
        dofRz = dofUx + 5;

        rows(end+1,1) = dofUx; cols(end+1,1) = colUx; vals(end+1,1) = 1.0;
        rows(end+1,1) = dofUx; cols(end+1,1) = colRy; vals(end+1,1) = -dz;
        rows(end+1,1) = dofUz; cols(end+1,1) = colUz; vals(end+1,1) = 1.0;
        rows(end+1,1) = dofUz; cols(end+1,1) = colRy; vals(end+1,1) = +dx;

        rows(end+1,1) = dofRy; cols(end+1,1) = colRy; vals(end+1,1) = 1.0;

        colUy = colIdx + 1; colIdx = colIdx + 1;  % UX-independent lateral offset
        rows(end+1,1) = dofUy; cols(end+1,1) = colUy; vals(end+1,1) = 1.0;

        colRx = colIdx + 1; colIdx = colIdx + 1;
        rows(end+1,1) = dofRx; cols(end+1,1) = colRx; vals(end+1,1) = 1.0;

        colRz = colIdx + 1; colIdx = colIdx + 1;
        rows(end+1,1) = dofRz; cols(end+1,1) = colRz; vals(end+1,1) = 1.0;
    end
end

% for floors below startFloor, keep all dofs uncoupled
if startFloor > 1
    for iFloor = 1:(startFloor-1)
        floorNodes = find(abs(nodeCoordinates(:,2)-levels(iFloor)) < 1e-9);
        for j = 1:numel(floorNodes)
            nodeId = floorNodes(j);
            for local = 1:6
                dofId = 6*(nodeId-1)+local;
                colIdx = colIdx + 1;
                rows(end+1,1) = dofId;
                cols(end+1,1) = colIdx;
                vals(end+1,1) = 1.0;
            end
        end
    end
end

% build reduced transform (full = T * reduced)
if isempty(rows)
    T = sparse([],[],[],GDof,0);
    dofToReduced = zeros(GDof,1);
else
    T = sparse(rows, cols, vals, GDof, colIdx);
    dofToReduced = zeros(GDof,1);
    for i = 1:GDof
        nz = find(T(i,:));
        if isempty(nz)
            error("applyRigidDiaphragmConstraints: unmapped DOF %d", i);
        end
        dofToReduced(i) = nz(1);
    end
end
stiffness = T' * stiffness * T;
mass = T' * mass * T;
force = T' * force;
end
