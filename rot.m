function R = rot(X, L, D)

% R = rot(X, L, D) returns a rotation matrix R (3x3) given
% the angle of rotation X in degree(radian) and the axis L.
% L must either be 'x', 'y', 'z' or a 3x1 vector.
% The rotation follows the right-hand screw rule.
% If D = 1, X is in degree; otherwise X is in radian.
% Default L = 'x', D = 1.
%
% see tx

% written by Kenneth Wong.
% last modified on 3-Nov-1999.

switch nargin
  case 1,
    L = 'x';
    D = 1;
  case 2,
    D = 1;
end

if D == 1
  X = X*pi/180;
end

ct = cos(X);
st = sin(X);

if ischar(L)
  switch L
    case 'x',
      R = [ 1   0   0
            0   ct -st
            0   st  ct];
    case 'y',
      R = [ ct  0   st
            0   1   0
           -st  0   ct];
    case 'z',
      R = [ ct -st  0
            st  ct  0
            0    0  1 ];
    otherwise,
      disp('L must be either x,y,z or a 3x1 vector!');
      R = zeros(3,3);
  end
  return
end


[m,n] = size(L);
if (m==3 & n==1)
  L = L/norm(L);
elseif (m==1 & n==3)
  L = L'/norm(L);
else
  disp('L must be either x,y,z or a 3x1 vector!');
  R = zeros(3,3);
  return
end

R = ct*eye(3)+st*Tx(L)+(1-ct)*L*L';

return
