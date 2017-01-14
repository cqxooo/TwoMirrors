function cost = getvz1D_new(data, vz)
fp = data.fp;
f = data.f;
cp = data.cp;
theta = data.theta;
lh = data.lh;
% if (lh(2)<0)
%     lhOrien = atan2(lh(1),-lh(2)); %ATAN2(Y,X);
% else
%     lhOrien = atan2(-lh(1),lh(2));
% end
u0 = distPnts(fp,vz(1:2));
K = [f u0;0 1];
i=cp(2,1);
j=cp(2,2);

vx(2,1) = (vz(2)*(i+j)-2*i*j)/(2*vz(2)-i-j);
vx(1,1) = -(lh(3)+lh(2)*vx(2))/lh(1);
vx(3,1) = 1;
kappa = real((cp(3,1)*vx(1)-cp(1,1)*vx(3))/(sqrt(-1)*(cp(3,1)*vz(1)-cp(1,1)*vz(3))));%
e = vx + kappa*tan(theta)*vz;
ev1 = e./e(3);
e = vx - kappa*tan(theta)*vz;
ev2 = e./e(3);
e = vx + kappa*tan(pi-2*theta)*vz;
ev21 = e./e(3);
e = vx - kappa*tan(pi-2*theta)*vz;
ev12 = e./e(3);

ev1D = distPnts(ev1, vz);
ev2D = distPnts(ev2, vz);
ev12D = distPnts(ev12, vz);
ev21D = distPnts(ev21, vz);
if ev1(1) < vz(1)
    ev1D = -ev1D;
end
if ev2(1) < vz(1)
    ev2D = -ev2D;
end
if ev12(1) < vz(1)
    ev12D = -ev12D;
end
if ev21(1) < vz(1)
    ev21D = -ev21D;
end

cost = 0;
%%%%%ev1 and ev2%%%%%%
H = K * [cos(-theta) sin(-theta);-sin(-theta) cos(-theta)]*inv(K);
ev2H = H*[ev2D;1];
ev2H = ev2H/ev2H(2);
cost = cost+abs(ev2H(1)-ev1D);
H = K * [cos(theta) sin(theta);-sin(theta) cos(theta)]*inv(K);
ev1H = H*[ev1D;1];
ev1H = ev1H/ev1H(2);
cost = cost+abs(ev1H(1)-ev2D);
%%%%%ev12 and ev21%%%%%%
H = K * [cos(-pi+2*theta) sin(-pi+2*theta);-sin(-pi+2*theta) cos(-pi+2*theta)]*inv(K);
ev12H = H*[ev12D;1];
ev12H = ev12H/ev12H(2);
cost = cost+abs(ev12H(1)-ev21D);

H = K * [cos(pi-2*theta) sin(pi-2*theta);-sin(pi-2*theta) cos(pi-2*theta)]*inv(K);
ev21H = H*[ev21D;1];
ev21H = ev21H/ev21H(2);
cost = cost+abs(ev21H(1)-ev12D);

EvalPrint('cost')
