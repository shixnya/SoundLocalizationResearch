function [basisset, pseudo_inverse, orthbasis] = cosyneBumpBasis(a, c, time, order, incl_const)
% unit of time is ms

if nargin < 5
    incl_const = 1;
end
if nargin < 4
    order = 12;
end
if nargin < 3
    time = 0:0.05:500;
end
if nargin < 2
    c = -5;
end    
if nargin < 2
    a = pi;
end

cs = 0;
basisset = zeros(length(time), order);
for i = 1:order
    phij = pi / 2 * (i - 1); % shifting term
    phase = a * log(time + c) - phij;
    phase(time + c <= 0) = pi;
    phase(phase < -pi) = pi;
    phase(phase > pi) = pi;
    b = 1/2 * cos(phase) + 1/2;
    cs = cs + b;
    basisset(:,i) = b;
end

if incl_const
    basisset = [ones(size(basisset, 1), 1), basisset];
end
pseudo_inverse = pinv(basisset);
orthbasis = orth(basisset);

