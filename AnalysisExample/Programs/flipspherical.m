function [ntheta, nphi] = flipspherical(theta, phi, mode)
% change the spherical coordinates.
if nargin < 3
    mode = 1; % front-z to top-z
end
x = sin(theta) .* cos(phi);
y = sin(theta) .* sin(phi);
z = cos(theta);

if mode == 1 % front-z to top-z
    nx = z;
    ny = x;
    nz = y;
else
    nx = y;
    ny = z;
    nz = x;
end
[nphi, ntheta] = cart2sph(nx, ny, nz);
ntheta = pi/2 - ntheta; % to convert matlab elev to physics theta.