% Total Least Squares
%   TLS(Data)
%   
%Return: [a,b,c] - line in a * x + b * y + c = 0 form
%                   where a ^ 2 + b ^ 2 = 1
%
%   TLS(X,Y) - approximates ALL points in array by one line

function line = TLS(Data)

  if any( size(Data) == 0)
      Line = [0, 0, 0];
      return;
  end

  X = Data(1, :);
  Y = Data(2, :);
  
  len = length(X);

  if size(X) ~= size(Y)
      X = X';
  end
    
  M = [ mid(X .^ 2) - mid(X) ^ 2, sum(X .* Y) / len - mid(X) * mid(Y);...
      sum(X .* Y) / len - mid(X) * mid(Y), mid(Y .^ 2) - mid(Y) ^ 2];
  [ev,tmp] = eig(M);
         
  ind = 2;
  if ErrFunc(X, Y, ev(:, 1)) < ErrFunc(X, Y, ev(:, 2))
      ind = 1;
  end
                  
  line = [ev(1,ind), ev(2,ind),...
         -ev(1,ind) * mid(X) - ev(2,ind) * mid(Y)];

return;

% Help function, calculates an error
function e = ErrFunc(X,Y,L)
maxE = 0;
e = 0;
c = -L(1) * mid(X) - L(2) * mid(Y);
for i = 1:length(X)
    e = e + ( L(1) * X(i) + L(2) * Y(i) + c )^2;
    if (L(1) * X(i) + L(2) * Y(i) + c )^2 > maxE
        maxE = (L(1) * X(i) + L(2) * Y(i) + c )^2;
    end;
end

return;

% Middle value of vector X
function l = mid(X)

if length(X) > 0
    l = sum(X) / length(X);
else 
    l = 0;
end;

return;