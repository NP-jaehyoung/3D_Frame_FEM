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

rows = [];
cols = [];
vals = [];
colIdx = 0;

% Floors below startFloor remain uncoupled (identity mapping)
if startFloor > 1
    for iFloor = 1:(startFloor-1)
        floorNodes = find(abs(nodeCoordinates(:,2)-levels(iFloor)) < 1e-9);
        for j = 1:numel(floorNodes)
            nodeId = floorNodes(j);
            for local = 1:6
                dofId = 6*(nodeId-1) + local;
                colIdx = colIdx + 1;
                rows(end+1,1) = dofId;
                cols(end+1,1) = colIdx;
                vals(end+1,1) = 1.0;
            end
        end
    end
end

% Floors from startFloor: rigid diaphragm in x-z plane + vertical offset handling
% Linearized rigid-body relation at node i (relative to floor centroid c):
% u_i = u_c + theta x r_i
% ux_i = ucx + (z-zc)*Ry - (y-yc)*Rz
% uy_i = ucy + (x-xc)*Rz - (z-zc)*Rx
% uz_i = ucz + (y-yc)*Rx - (x-xc)*Ry
for iFloor = startFloor:numFloors
    floorNodes = find(abs(nodeCoordinates(:,2)-levels(iFloor)) < 1e-9);
    if isempty(floorNodes)
        continue;
    end

    cx = mean(nodeCoordinates(floorNodes,1));
    cy = mean(nodeCoordinates(floorNodes,2));
    cz = mean(nodeCoordinates(floorNodes,3));

    colUx = colIdx + 1; colIdx = colIdx + 1;
    colUy = colIdx + 1; colIdx = colIdx + 1;
    colUz = colIdx + 1; colIdx = colIdx + 1;
    colRx = colIdx + 1; colIdx = colIdx + 1;
    colRy = colIdx + 1; colIdx = colIdx + 1;
    colRz = colIdx + 1; colIdx = colIdx + 1;

    for j = 1:numel(floorNodes)
        nodeId = floorNodes(j);
        dx = nodeCoordinates(nodeId,1) - cx;
        dy = nodeCoordinates(nodeId,2) - cy;
        dz = nodeCoordinates(nodeId,3) - cz;

        dofUx = 6*(nodeId-1) + 1;
        dofUy = dofUx + 1;
        dofUz = dofUx + 2;
        dofRx = dofUx + 3;
        dofRy = dofUx + 4;
        dofRz = dofUx + 5;

        % ux_i = ucx + dz*Ry - dy*Rz
        rows(end+1,1) = dofUx; cols(end+1,1) = colUx; vals(end+1,1) = 1.0;
        rows(end+1,1) = dofUx; cols(end+1,1) = colRy; vals(end+1,1) =  dz;
        rows(end+1,1) = dofUx; cols(end+1,1) = colRz; vals(end+1,1) = -dy;

        % uy_i = ucy + dx*Rz - dz*Rx
        rows(end+1,1) = dofUy; cols(end+1,1) = colUy; vals(end+1,1) = 1.0;
        rows(end+1,1) = dofUy; cols(end+1,1) = colRz; vals(end+1,1) =  dx;
        rows(end+1,1) = dofUy; cols(end+1,1) = colRx; vals(end+1,1) = -dz;

        % uz_i = ucz + dy*Rx - dx*Ry
        rows(end+1,1) = dofUz; cols(end+1,1) = colUz; vals(end+1,1) = 1.0;
        rows(end+1,1) = dofUz; cols(end+1,1) = colRx; vals(end+1,1) =  dy;
        rows(end+1,1) = dofUz; cols(end+1,1) = colRy; vals(end+1,1) = -dx;

        % rigid rotations
        rows(end+1,1) = dofRx; cols(end+1,1) = colRx; vals(end+1,1) = 1.0;
        rows(end+1,1) = dofRy; cols(end+1,1) = colRy; vals(end+1,1) = 1.0;
        rows(end+1,1) = dofRz; cols(end+1,1) = colRz; vals(end+1,1) = 1.0;
    end
end

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
