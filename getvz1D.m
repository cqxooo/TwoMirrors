function cost = getvz1D(data, vz)
fp = data.fp;
cp = data.cp;
theta = data.theta;
lh = data.lh;
if (lh(2)<0)
    lhOrien = atan2(lh(1),-lh(2)); %ATAN2(Y,X);
else
    lhOrien = atan2(-lh(1),lh(2));
end
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
i = cp(:,1);
j = cp(:,2);
a = (i(1)+j(1)-2*vz(1))/(2*cos(lhOrien));
b = (i(1)-j(1))/(2*sqrt(-1)*cos(lhOrien));

% a = (i(2)+j(2)-2*vz(2))/(2*sin(lhOrien));
% b = (i(2)-j(2))/(2*sqrt(-1)*sin(lhOrien));
c = [a+sqrt(-1)*b;1];
V = [c conj(c)];

cost = 0;
%%%%%ev1 and ev2%%%%%%
d = exp(-sqrt(-1)*theta);
D = [d 0;0 conj(d)];
H = real(V*D*inv(V));
H = H/H(2,1);

ev2H = H*[ev2D;1];
ev2H = ev2H/ev2H(2);
cost = cost+abs(ev2H(1)-ev1D);
% ev1H = H*[ev1D;1];
% ev1H = ev1H/ev1H(2);
% cost = cost+abs(ev1H(1)-ev2D);
%%%%%ev12 and ev21%%%%%%
d = exp(-sqrt(-1)*(pi-2*theta));
D = [d 0;0 conj(d)];
H = real(V*D*inv(V));
H = H/H(2,1);

% ev12H = H*[ev12D;1];
% ev12H = ev12H/ev12H(2);
% cost = cost+abs(ev12H(1)-ev21D);
ev21H = H*[ev21D;1];
ev21H = ev21H/ev21H(2);
cost = cost+abs(ev21H(1)-ev12D);
EvalPrint('cost')
