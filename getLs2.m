function [costSum ls] = getLs2(polyBoundaryVecCell, cp, theta, lh, vz, T)
pbv=polyBoundaryVecCell{1};
s = size(pbv,2);
O =  T*[pbv; ones(1,s)];
Or = [O(1,:)./O(3,:); O(2,:)./O(3,:); O(3,:)./O(3,:)];
pbv=polyBoundaryVecCell{2};
s = size(pbv,2);
O =  T*[pbv; ones(1,s)];
Ov1 = [O(1,:)./O(3,:); O(2,:)./O(3,:); O(3,:)./O(3,:)];
pbv=polyBoundaryVecCell{3};
s = size(pbv,2);
O =  T*[pbv; ones(1,s)];
Ov2 = [O(1,:)./O(3,:); O(2,:)./O(3,:); O(3,:)./O(3,:)];
pbv=polyBoundaryVecCell{4};
s = size(pbv,2);
O =  T*[pbv; ones(1,s)];
Ov12 = [O(1,:)./O(3,:); O(2,:)./O(3,:); O(3,:)./O(3,:)];
pbv=polyBoundaryVecCell{5};
s = size(pbv,2);
O =  T*[pbv; ones(1,s)];
Ov21 = [O(1,:)./O(3,:); O(2,:)./O(3,:); O(3,:)./O(3,:)];

vz=vz./vz(3);
i=cp(2,1);
j=cp(2,2);
v = vz(2);
l = [lh(2); -lh(1); -lh(2)*vz(1)+lh(1)*vz(2)];
ls = -sign(l(3))*l/norm(l(1:2));


vx(2) = (v*(i+j)-2*i*j)/(2*v-i-j);
vx(1) = -(lh(3)+lh(2)*vx(2))/lh(1);
vx(3) = 1;
vx = vx';

cost = [];

c = Tx(ls)*lh;
kappa = real((cp(3,1)*vx(1)-cp(1,1)*vx(3))/(sqrt(-1)*(cp(1,1)*c(3)-cp(3,1)*c(1))));%
e = vx - kappa*tan(theta)*c;
ev1 = e./e(3);
e = vx + kappa*tan(theta)*c;
ev2 = e./e(3);
e = vx - kappa*tan(pi-2*theta)*Tx(ls)*lh;
ev21 = e./e(3);
e = vx + kappa*tan(pi-2*theta)*Tx(ls)*lh;
ev12 = e./e(3);

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Ov12',ev12);
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(Ov21',ev21);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
pointUp = cross(linelUp,linerUp);
pointUp = pointUp./pointUp(3);
pointBot = cross(linelBot,linerBot);
pointBot = pointBot./pointBot(3);
t1 = abs(pointUp'*ls);
t2 = abs(pointBot'*ls);
cost=[cost t1 t2];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Ov12',ev1);
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(Or',ev2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
pointUp = cross(linelUp,linerUp);
pointUp = pointUp./pointUp(3);
pointBot = cross(linelBot,linerBot);
pointBot = pointBot./pointBot(3);
t1 = abs(pointUp'*ls);
t2 = abs(pointBot'*ls);
cost=[cost t1 t2];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Or',ev1);
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(Ov21',ev2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
pointUp = cross(linelUp,linerUp);
pointUp = pointUp./pointUp(3);
pointBot = cross(linelBot,linerBot);
pointBot = pointBot./pointBot(3);
t1 = abs(pointUp'*ls);
t2 = abs(pointBot'*ls);
cost=[cost t1 t2];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Ov1',ev1);
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(Ov2',ev2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
pointUp = cross(linelUp,linerUp);
pointUp = pointUp./pointUp(3);
pointBot = cross(linelBot,linerBot);
pointBot = pointBot./pointBot(3);
t1 = abs(pointUp'*ls);
t2 = abs(pointBot'*ls);
cost=[cost t1 t2];

costSum = sum(cost);
EvalPrint('costSum')