function d = distPnts(a,b)
% obtain the distance between 2 points
if (size(a,1)==2)
    d = sqrt((a(1)-b(1))*(a(1)-b(1)) + (a(2)-b(2))*(a(2)-b(2)));
elseif (size(a,1)==3)
    d = sqrt((a(1)-b(1))*(a(1)-b(1)) + (a(2)-b(2))*(a(2)-b(2))+(a(3)-b(3))*(a(3)-b(3)));
end
    
return;