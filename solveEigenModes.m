function [modes,eigenvalues] = solveEigenModes(GDof,prescribedDof, ...
    stiffness,mass,maxEigenvalues)
activeDof = setdiff((1:GDof)', prescribedDof);
if maxEigenvalues == 0
    [V,D] = eig(stiffness(activeDof,activeDof),...
        mass(activeDof,activeDof));
else
    [V,D] = eigs(stiffness(activeDof,activeDof),...
        mass(activeDof,activeDof),maxEigenvalues,'smallestabs');
end
eigenvalues = diag(D);
modes = zeros(GDof,length(eigenvalues));
modes(activeDof,:) = V;
end
