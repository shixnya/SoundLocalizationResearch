function arrowaxis(startp, endp, varargin)
% place text in a relative position of the axis
axp = get(gca, 'Position');
startp
endp
axp
x = [axp(1) + startp(1) * axp(3), axp(1) + endp(1) * axp(3)];
y = [axp(2) + startp(2) * axp(4), axp(2) + endp(2) * axp(4)];
x
y
annotation('arrow', x, y, varargin{1:end});
