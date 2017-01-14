function angle = getTheta(cp,v)
v1 = v(:,1);
v2 = v(:,2);
v3 = v(:,3);
v4 = v(:,4);
theta(1) = log(crossRatio(v1, v2, cp(:,1), cp(:,2)))/2*sqrt(-1);
theta(2) = log(crossRatio(v3, v1, cp(:,1), cp(:,2)))/2*sqrt(-1);
theta(3) = log(crossRatio(v2, v4, cp(:,1), cp(:,2)))/2*sqrt(-1);
angle = mean(real(theta));