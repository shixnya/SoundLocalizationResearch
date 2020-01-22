function colors = bluered()

colors = ones(256, 3);

grad = (1:128) / 128;

colors(1:128, 1) = grad;
colors(1:128, 2) = grad;

colors(256:-1:129, 3) = grad;
colors(256:-1:129, 2) = grad;
