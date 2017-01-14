function D = kenDist(A, B)

% D = kenDist(A, B) returns the distances D between
% each points in A and B.
%
% both A and B are nx2/nx3 matrices containing the
% (x,y)/(x,y,z) coordinates of 2 sets of n points.
%
% or
%
% A is a nx2/nx3 matrix and B is a 1x2/1x3 vector,
% or vice versa, containing the (x,y)/(x,y,z)
% coordinates of points. In this case, D contains
% a list of distances between the single point
% and each point in the list.
%

D = [];

[m1,n1] = size(A);
[m2,n2] = size(B);

if (n1 ~= n2 | (n1 ~=2 & n1 ~=3))
  disp('Input matrices must be nx2/nx3!')
  return
end

if (m1 > 1 & m2 == 1) 
  B = repmat(B, m1, 1);
elseif (m1 == 1 & m2 > 1)
  A = repmat(A, m2, 1);
else
  s = min([m1 m2]);
  A = A(1:s,:);
  B = B(1:s,:);
end

d = A - B;

[m,n] = size(d);

switch n
case 2,
  d = d(:,1).^2 + d(:,2).^2;
case 3,
  d = d(:,1).^2 + d(:,2).^2 + d(:,3).^2;
end

D = sqrt(d);

return