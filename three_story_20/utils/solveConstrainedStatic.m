function displacements = solveConstrainedStatic(GDof, prescribedDof, stiffness, force)
activeDof = setdiff((1:GDof)', [prescribedDof]);
U = stiffness(activeDof,activeDof)\force(activeDof);
displacements = zeros(GDof,1);
displacements(activeDof) = U;
end
