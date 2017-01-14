function ls = Initvzls_new(polyBoundaryVecCell, cp, theta, lh, v, T)
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

v = [ v(1,:)./v(3,:); v(2,:)./v(3,:); v(3,:)./v(3,:)];
a = distPnts(v(:,1),v(:,4));
b = distPnts(v(:,3),v(:,4));
c = distPnts(v(:,2),v(:,3));
lhDir = (v(1:2,2)- v(1:2,1))/norm(v(1:2,2)- v(1:2,1));
fp = v(1:2,1) + 1/2*(2*a+2*b+c)*a*(a+b+c)/(a^2+a*b+c^2+c*b+a*c) * lhDir;
fp = [fp;1];
i=cp(2,1);
j=cp(2,2);
vx(2,1) = (fp(2)*(i+j)-2*i*j)/(2*fp(2)-i-j);
vx(1,1) = -(lh(3)+lh(2)*vx(2))/lh(1);
vx(3,1) = 1;
kappa = real((cp(3,1)*vx(1)-cp(1,1)*vx(3))/(sqrt(-1)*(cp(3,1)*fp(1)-cp(1,1)*fp(3))));%
e = vx + kappa*tan(theta)*fp;
ev1 = e./e(3);
e = vx - kappa*tan(theta)*fp;
ev2 = e./e(3);
e = vx + kappa*tan(pi-2*theta)*fp;
ev21 = e./e(3);
e = vx - kappa*tan(pi-2*theta)*fp;
ev12 = e./e(3);

p=fp;
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
p1=cross(linelUp, linerUp);
p2=cross(linelBot,linerBot);
p1 = p1/p1(3);
p2 = p2/p2(3);
p=[p p1 p2];

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
p1=cross(linelUp, linerUp);
p2=cross(linelBot,linerBot);
p1 = p1/p1(3);
p2 = p2/p2(3);
p=[p p1 p2];

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
p1=cross(linelUp, linerUp);
p2=cross(linelBot,linerBot);
p1 = p1/p1(3);
p2 = p2/p2(3);
p=[p p1 p2];

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
p1=cross(linelUp, linerUp);
p2=cross(linelBot,linerBot);
p1 = p1/p1(3);
p2 = p2/p2(3);
p=[p p1 p2];

[u d v] = svd(p');
v=v(:,3);
ls= -sign(v(3))*v/norm(v(1:2));

fig=figure(5);
ShowPoly( Or, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
hold on
axis ij;
axis equal;
ShowPoly( Ov1, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15); 
ShowPoly( Ov2, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
ShowPoly( Ov12, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
ShowPoly( Ov21, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15);  
drawLine(ls, 'g-');
close(fig)


