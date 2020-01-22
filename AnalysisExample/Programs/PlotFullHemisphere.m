function h = PlotFullHemisphere(counts)
% plot the 17 x 5 hemispheric pattern in the spherecal plot.



% first convert the format to 9 x 10 space.
newcounts = [counts(:, 9:end); counts(end:-1:1, 9:-1:1)];

[oazim, oelev] = meshgrid((0:18:144) * pi / 180, (0:20:180) * pi / 180);

x = cos(oazim);
y = -sin(oazim) .* cos(oelev);
z = sin(oazim) .* sin(oelev);

h = surf(x, y, z, newcounts);
shading interp
axis image
campos([1 -1 1])