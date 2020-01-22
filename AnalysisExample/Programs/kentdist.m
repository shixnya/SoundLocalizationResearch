function kent = kentdist(p, x)
% provide unnormalized Kent distribution.
% p is 5 parameters for the kent distribution
% modified on 6/10/19 I should divide the peak by e(kappa)...

sph2cart = @(theta, phi) cat(3, sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta));

kappa = p(1);
beta = p(2);
theta = p(3);
phi = p(4);
alpha = p(5);
amp = p(6);
baseline = p(7);

theta_x = x{1};
phi_x = x{2};

units = sphericalUnit(theta, phi);
gamma1 = units(:, 1);
gamma2 = rodrot(units(:, 2), units(:, 1), alpha);
gamma3 = rodrot(units(:, 3), units(:, 1), alpha);

gamma1 = permute(gamma1, [2, 3, 1]);
gamma2 = permute(gamma2, [2, 3, 1]);
gamma3 = permute(gamma3, [2, 3, 1]);

gamma1 = repmat(gamma1, [size(theta_x), 1]);
gamma2 = repmat(gamma2, [size(theta_x), 1]);
gamma3 = repmat(gamma3, [size(theta_x), 1]);

xyz = sph2cart(theta_x, phi_x);
kent = amp * exp(-kappa) * exp(kappa * dot(gamma1, xyz, 3) + ...
    beta * kappa * (dot(gamma2, xyz, 3).^2 - dot(gamma3, xyz, 3).^2)) + baseline;


%figure(15112);
%surf(x, y, z, kent);
%shading interp
%axis image
