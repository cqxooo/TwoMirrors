function C = fitBSpline(P, N)

% C = fitBSpline(P, N) returns a list of control points
% C (?x2) for a closed splines fitted to the points P (?x2).
% each segment is fitted to N points in the list.

len = size(P,1);
Nsp = ceil(len/N);

A = zeros(len,Nsp);

fixed_matrix = [-1  3 -3 1
                 3 -6  3 0
                -3  0  3 0
                 1  4  1 0]/6;

samples = [];
samp_index = [0:N-1]/N;

for s = samp_index
  samples = [samples; [s^3 s^2 s 1]*fixed_matrix];
end

% construct the samples for the first Nsp-1 splines
for j=1:Nsp-1
  A((j-1)*N + (1:N), cycle(j:j+3,Nsp)) = samples;
end

Nl = mod(len,N);
if Nl == 0
  Nl = N;
end

% construct the samples for the last spline
samples = [];
samp_index = [0:Nl-1]/Nl;
for s = samp_index
  samples = [samples; [s^3 s^2 s 1]*fixed_matrix];
end
A((len-Nl+1):len, [Nsp 1 2 3]) = samples;

pinvA = pinv(A);

C = pinvA*P;

return