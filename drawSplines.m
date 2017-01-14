function [H, X] = drawSplines(B, N, S, C)

% [H, X] = drawSplines(B, N, S, C) draws the B-Splines with the given 
% control points B(?,2) with N samples per segments with style S.
% If C = 1, the B-spline will be closed, otherwise it is open.
% Default N = 10, S = 'r-' and C = 0.
%
% The sample points are returned in X ((ns*N)x3) and the handle
% of the B-spline is returned in H.
%

switch nargin
case 1,
  N = 10;
  S = 'r-';
  C = 0;
case 2,
  S = 'r-';
  C = 0;
case 3,
  C = 0;
end

fixed_matrix = [-1  3 -3  1
                 3 -6  3  0
                -3  0  3  0
                 1  4  1  0]/6;

Ncp = size(B,1);

if C == 1
  Ns = Ncp;
else
  Ns = Ncp - 3;
end

A = zeros(Ns*N,Ncp);

samples = [];
samp_index = [0:(N-1)]/N;
for s = samp_index
  samples = [samples; [s^3 s^2 s 1]*fixed_matrix];
end

for j=1:Ns
  A((j-1)*N + (1:N), cycle(j:j+3, Ncp)) = samples;
end

X = [A*B ones(Ns*N,1)];

if C==1
    H = plot([X(:,1); X(1,1)],[X(:,2); X(1,2)],S,'LineWidth',2);
else
  X = [X; [[0 1 4 1]*B((Ncp-3:Ncp),:)/6 1]];
  H = plot(X(:,1),X(:,2),S);
end

return



