function C = cycle(I,N)

% C = cycle(I,N) returns the number C within
% the range 1 to N given I by cycling through 
% the range. 
%

C = mod(I-1,N) + 1;

return