function x_ = getPerPtLine(x, l)
  
% x_ = getPerPtLine(x, l) returns the orthogonal projection of
% the point x on the line l
%
  
l = -sign(l(3))*l/norm(l(1:2));
d = l(1)*x(1)+l(2)*x(2)+l(3);
x_ = x - d*l(1:2);

return
