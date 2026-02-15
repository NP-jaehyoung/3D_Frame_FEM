% ................................................................
function outputDisplacementsReactions...
    (displacements,stiffness,GDof,prescribedDof,force)
% output of displacements and reactions in tabular form
if nargin < 5 || isempty(force)
    force = zeros(GDof,1);
end

% GDof: total number of degrees of freedom of the problem

% displacements
disp('Displacements')
jj = 1:GDof; format
[jj' displacements]

% reactions
F = stiffness*displacements;
reactions = F(prescribedDof) - force(prescribedDof);
disp('reactions')
[prescribedDof reactions]

end
