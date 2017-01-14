function ls = LineRANSAC(points,num)
len = size(points,2);


Idx = len * rand(2,num);
Idx = ceil(Idx);
minMedian = 10000;
for i=1:num
    l = cross(points(:,Idx(1,i)),points(:,Idx(2,i)));
    l = -sign(l(3))*l/norm(l(1:2));

    d = abs(l'*points);
    tMedian = median(d);
    if tMedian < minMedian
        minMedian = tMedian;
        ls = l;
    end
end
return;