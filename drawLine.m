function H = drawLine(L, S)

% H = drawLine(L, S) draws the L with style S across
% the whole figure and the handle H is returned.
%

aux = axis;
if (abs(L(1)) > abs(L(2)))
  y0 = aux(3);
  x0 = (-L(2)*y0 - L(3))/L(1);

  yf = aux(4);
  xf = (-L(2)*yf - L(3))/L(1);
else
  x0 = aux(1);
  y0 = (-L(1)*x0 - L(3))/L(2);

  xf = aux(2);
  yf = (-L(1)*xf - L(3))/L(2);
end

H = plot([x0, xf], [y0, yf], S, 'LineWidth',2);
% H = plot([x0, xf], [y0, yf], S);
axis(aux);

return

