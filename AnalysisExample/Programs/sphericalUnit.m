function unitvecs = sphericalUnit(theta, phi)
% this function gives a unit vectors of spherical coordinates.
% the notation is based on Arfken
% theta is polar angle
% phi is azimuthal angle.

st = sin(theta);
ct = cos(theta);
sp = sin(phi);
cp = cos(phi);

unitvecs = [st * cp,  ct * cp, -sp;...
            st * sp,  ct * sp,  cp;...
            ct,      -st,      0];


    
