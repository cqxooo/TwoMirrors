function [lh,v0, v1, residual]=getRefinedLh(polyBoundaryVecCell,lh,v,T)
% this function tries to refine the coeficent ktan(theta/2).
% data structure has ls, vx,lh,Cl, Cr.
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
v0 = v;

theta0 = atan2(lh(2),lh(1));
d0 = -lh(3);
rho1 = atan2(v(2,1),v(1,1));
rho2 = atan2(v(2,2),v(1,2));
rho3 = atan2(v(2,3),v(1,3));
rho4 = atan2(v(2,4),v(1,4));
param0 = [theta0 d0 rho1 rho2 rho3 rho4];
data = [];
data.Or = Or;
data.Ov1 = Ov1;
data.Ov2 = Ov2;
data.Ov12 = Ov12;
data.Ov21 = Ov21;

cost=costFunctionDistanceSilhouettes(param0,data);

% set the options for lsqnonlin
options = optimset('lsqnonlin');
options = optimset(options, 'display', 'iter');
options = optimset(options,'TolFun',1e-20);
options = optimset(options,'TolX',1e-20);
options = optimset(options,'MaxIter',30);
options = optimset(options,'MaxFunEvals',800);
% set the options for lsqnonlin
[param,resnorm,residual] = lsqnonlin(@costFunctionDistanceSilhouettes, param0,...
              [],[],options,data);
theta = param(1);
d = param(2);
rho1 = param(3);
rho2 = param(4);
rho3 = param(5);
rho4 = param(6);

l = [cos(theta); sin(theta); -d];
lh = -sign(l(3))*l/norm(l(1:2));
r = rho1;
v = [cos(r); sin(r); -(l(1)*cos(r)+l(2)*sin(r))/l(3)];
v1 = v/v(3);
r = rho2;
v = [cos(r); sin(r); -(l(1)*cos(r)+l(2)*sin(r))/l(3)];
v2 = v/v(3);
r = rho3;
v = [cos(r); sin(r); -(l(1)*cos(r)+l(2)*sin(r))/l(3)];
v3 = v/v(3);
r = rho4;
v = [cos(r); sin(r); -(l(1)*cos(r)+l(2)*sin(r))/l(3)];
v4 = v/v(3);
v = [v1 v2 v3 v4];
v1 = v;

fig=figure(3);
ShowPoly( Or, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
hold on
axis ij;
axis equal;   
ShowPoly( Ov1, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15);  
ShowPoly( Ov2, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
ShowPoly( Ov12, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
ShowPoly( Ov21, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15);  
plot(v1(1),v1(2),'k+');
plot(v2(1),v2(2),'k+');
plot(v3(1),v3(2),'k+');
plot(v4(1),v4(2),'k+');
drawLine(lh,'g-');
close(fig);
return;


function cost=costFunctionDistanceSilhouettes(param,data)
theta = param(1);
d = param(2);
rho1 = param(3);
rho2 = param(4);
rho3 = param(5);
rho4 = param(6);

l = [cos(theta); sin(theta); -d];
lh = -sign(l(3))*l/norm(l(1:2));
r = rho1;
v = [cos(r); sin(r); -(lh(1)*cos(r)+lh(2)*sin(r))/lh(3)];
v1 = v/v(3);
r = rho2;
v = [cos(r); sin(r); -(lh(1)*cos(r)+lh(2)*sin(r))/lh(3)];
v2 = v/v(3);
r = rho3;
v = [cos(r); sin(r); -(lh(1)*cos(r)+lh(2)*sin(r))/lh(3)];
v3 = v/v(3);
r = rho4;
v = [cos(r); sin(r); -(lh(1)*cos(r)+lh(2)*sin(r))/lh(3)];
v4 = v/v(3);
v = [v1 v2 v3 v4];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Or',v1);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(data.Ov1',v1);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
t1 = abs(pointrUp'*linelUp)   ;
t2 = abs(pointrBot'*linelBot) ;
t3 = abs(pointlUp'*linerUp)   ;
t4 = abs(pointlBot'*linerBot);
cost = [t1 t2 t3 t4];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Ov2',v1);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(data.Ov12',v1);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
t1 = abs(pointrUp'*linelUp)   ;
t2 = abs(pointrBot'*linelBot) ;
t3 = abs(pointlUp'*linerUp)   ;
t4 = abs(pointlBot'*linerBot);
cost = [cost t1 t2 t3 t4];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Or',v2);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(data.Ov2',v2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
t1 = abs(pointrUp'*linelUp)   ;
t2 = abs(pointrBot'*linelBot) ;
t3 = abs(pointlUp'*linerUp)   ;
t4 = abs(pointlBot'*linerBot);
cost = [cost t1 t2 t3 t4];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Ov1',v2);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(data.Ov21',v2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
t1 = abs(pointrUp'*linelUp)   ;
t2 = abs(pointrBot'*linelBot) ;
t3 = abs(pointlUp'*linerUp)   ;
t4 = abs(pointlBot'*linerBot);
cost = [cost t1 t2 t3 t4];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Ov1',v3);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(data.Ov12',v3);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
t1 = abs(pointrUp'*linelUp)   ;
t2 = abs(pointrBot'*linelBot) ;
t3 = abs(pointlUp'*linerUp)   ;
t4 = abs(pointlBot'*linerBot);
cost = [cost t1 t2 t3 t4];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Ov2',v4);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(data.Ov21',v4);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
t1 = abs(pointrUp'*linelUp)   ;
t2 = abs(pointrBot'*linelBot) ;
t3 = abs(pointlUp'*linerUp)   ;
t4 = abs(pointlBot'*linerBot);
cost = [cost t1 t2 t3 t4];
return;