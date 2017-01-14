function cr = crossRatio(u, v, x, y)
% obtain the cross ratio between 4 points

cr = (distPnts(u,x)/distPnts(u,y)) * (distPnts(v,y)/distPnts(v,x));

return;