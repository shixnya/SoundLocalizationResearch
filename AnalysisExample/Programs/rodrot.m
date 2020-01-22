function rotated = rodrot(targetvector, rotationaxis, angle)
% this function does rotation of a vector in 3d space accordingly to
% Rodrigues rotation formula.

v = targetvector;
k = rotationaxis;
t = angle;

rotated = v * cos(t) + cross(k, v) * sin(t) + k * (k' * v) * (1 - cos(t));