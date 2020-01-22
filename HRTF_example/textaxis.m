function textaxis(sometext, relpos, varargin)
% place text in a relative position of the axis
axp = get(gca, 'Position');
resp = [axp(1) + relpos(1) * axp(3), axp(2), 1, relpos(2) * axp(4)];
annotation('textbox', resp, 'String', sometext, 'EdgeColor', 'none', varargin{1:end});
