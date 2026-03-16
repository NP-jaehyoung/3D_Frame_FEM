function [A, Iy, Iz, J, G] = frameSectionProperties(width, depth, E, nu)
% Compute rectangular 3D frame section properties used by the sample models.

if nargin < 4 || isempty(nu)
    nu = 0.167;
end

if width <= 0 || depth <= 0
    error("width and depth must be positive.");
end

A = width * depth;
Iy = width * depth^3 / 12;
Iz = depth * width^3 / 12;

if width >= depth
    ar = depth / width;
    J = depth * width^3 * (1/3 - 0.21*ar + 0.063*ar^5);
else
    ar = width / depth;
    J = width * depth^3 * (1/3 - 0.21*ar + 0.063*ar^5);
end

if J <= 0
    error("Invalid torsional constant J (non-positive). Check width/depth.");
end

G = E / (2 * (1 + nu));
end
