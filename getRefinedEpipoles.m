function data1 = getRefinedEpipoles(polyBoundaryVecCell, lh, vz, cp, theta, T)
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

data.Or = Or;
data.Ov1 = Ov1;
data.Ov2 = Ov2;
data.Ov12 = Ov12;
data.Ov21 = Ov21;
data.lh=lh;
data.cp=cp;
data.theta=theta;

phi = atan2(vz(2),vz(1));

param0 = phi;
% set the options for lsqnonlin
options = optimset('lsqnonlin');
options = optimset(options, 'display', 'iter');
options = optimset(options,'TolFun',1e-20);
options = optimset(options,'TolX',1e-20);
options = optimset(options,'MaxIter',80);
options = optimset(options,'MaxFunEvals',100);
% set the options for lsqnonlin
cost=costFunction_ls(param0,data);
[param,resnorm,residual] = lsqnonlin(@costFunction_ls, param0,...
    [],[],options,data);

phi = param(1);
vz = [cos(phi); sin(phi); -(lh(1)*cos(phi)+lh(2)*sin(phi))/lh(3)];
vz = vz/vz(3);
i=cp(2,1);
j=cp(2,2);
vx(2,1) = (vz(2)*(i+j)-2*i*j)/(2*vz(2)-i-j);
vx(1,1) = -(lh(3)+lh(2)*vx(2))/lh(1);
vx(3,1) = 1;

ls = getLs_cost(data, vz);

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
close(fig);

l = T'*ls;
ls = -sign(l(3))*l/norm(l(1:2));%
l = T'*lh;
lh = -sign(l(3))*l/norm(l(1:2));%
v = inv(T)*vx;
vx = v./v(3);%
c = cp;
v = inv(T)*c;
cp = [v(1,:)./v(3,:);v(2,:)./v(3,:);v(3,:)./v(3,:)] ;%

data1.ls = ls;
data1.lh = lh;
data1.vx = vx;
data1.cp = cp;
data1.Or=Or;
data1.Ov1=Ov1;
data1.Ov2=Ov2;
data1.Ov12=Ov12;
data1.Ov21=Ov21;

return;


function cost=costFunction_ls(param,data)

lh = data.lh;
theta = data.theta;
cp = data.cp;

phi = param(1);
vz = [cos(phi); sin(phi); -(lh(1)*cos(phi)+lh(2)*sin(phi))/lh(3)];
vz = vz/vz(3);

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

ls = getLs_cost(data, vz);

W = eye(3)-2*vx*ls'/(vx'*ls);
cost = [];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Ov1',ev1);% note the order is changed because of the position
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(data.Ov2',ev2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end

pointUp = cross(linelUp,linerUp);
pointUp = pointUp./pointUp(3);
pointBot = cross(linelBot,linerBot);
pointBot = pointBot./pointBot(3);
linelUpW = inv(W)'*linelUp;
linelUpW = -sign(linelUpW(3))*linelUpW/norm(linelUpW(1:2));
linelBotW = inv(W)'*linelBot;
linelBotW = -sign(linelBotW(3))*linelBotW/norm(linelBotW(1:2));
linerUpW = inv(W)'*linerUp;
linerUpW = -sign(linerUpW(3))*linerUpW/norm(linerUpW(1:2));
linerBotW = inv(W)'*linerBot;
linerBotW = -sign(linerBotW(3))*linerBotW/norm(linerBotW(1:2));
t1 = abs(pointlUp'*linerUpW)   ;
t2 = abs(pointlBot'*linerBotW) ;
t3 = abs(pointrUp'*linelUpW)   ;
t4 = abs(pointrBot'*linelBotW);
t5 = abs(pointUp'*ls);
t6 = abs(pointBot'*ls);
cost = [cost t1 t2 t3 t4 t5 t6];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Ov12',ev1);
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(data.Or',ev2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
pointUp = cross(linelUp,linerUp);
pointUp = pointUp./pointUp(3);
pointBot = cross(linelBot,linerBot);
pointBot = pointBot./pointBot(3);
linelUpW = inv(W)'*linelUp;
linelUpW = -sign(linelUpW(3))*linelUpW/norm(linelUpW(1:2));
linelBotW = inv(W)'*linelBot;
linelBotW = -sign(linelBotW(3))*linelBotW/norm(linelBotW(1:2));
linerUpW = inv(W)'*linerUp;
linerUpW = -sign(linerUpW(3))*linerUpW/norm(linerUpW(1:2));
linerBotW = inv(W)'*linerBot;
linerBotW = -sign(linerBotW(3))*linerBotW/norm(linerBotW(1:2));
t1 = abs(pointlUp'*linerUpW)   ;
t2 = abs(pointlBot'*linerBotW) ;
t3 = abs(pointrUp'*linelUpW)   ;
t4 = abs(pointrBot'*linelBotW);
t5 = abs(pointUp'*ls);
t6 = abs(pointBot'*ls);
cost = [cost t1 t2 t3 t4 t5 t6];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Or',ev1);
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(data.Ov21',ev2);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
pointUp = cross(linelUp,linerUp);
pointUp = pointUp./pointUp(3);
pointBot = cross(linelBot,linerBot);
pointBot = pointBot./pointBot(3);
linelUpW = inv(W)'*linelUp;
linelUpW = -sign(linelUpW(3))*linelUpW/norm(linelUpW(1:2));
linelBotW = inv(W)'*linelBot;
linelBotW = -sign(linelBotW(3))*linelBotW/norm(linelBotW(1:2));
linerUpW = inv(W)'*linerUp;
linerUpW = -sign(linerUpW(3))*linerUpW/norm(linerUpW(1:2));
linerBotW = inv(W)'*linerBot;
linerBotW = -sign(linerBotW(3))*linerBotW/norm(linerBotW(1:2));
t1 = abs(pointlUp'*linerUpW)   ;
t2 = abs(pointlBot'*linerBotW) ;
t3 = abs(pointrUp'*linelUpW)   ;
t4 = abs(pointrBot'*linelBotW);
t5 = abs(pointUp'*ls);
t6 = abs(pointBot'*ls);
cost = [cost t1 t2 t3 t4 t5 t6];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(data.Ov21',ev21);% note the order is changed because of the position
[linerBot,linerUp,pointrBot,pointrUp,flagr] = getOuterTengents(data.Ov12',ev12);
if flagl == 0 || flagr == 0
    cost = 1000;
    return;
end
pointUp = cross(linelUp,linerUp);
pointUp = pointUp./pointUp(3);
pointBot = cross(linelBot,linerBot);
pointBot = pointBot./pointBot(3);
linelUpW = inv(W)'*linelUp;
linelUpW = -sign(linelUpW(3))*linelUpW/norm(linelUpW(1:2));
linelBotW = inv(W)'*linelBot;
linelBotW = -sign(linelBotW(3))*linelBotW/norm(linelBotW(1:2));
linerUpW = inv(W)'*linerUp;
linerUpW = -sign(linerUpW(3))*linerUpW/norm(linerUpW(1:2));
linerBotW = inv(W)'*linerBot;
linerBotW = -sign(linerBotW(3))*linerBotW/norm(linerBotW(1:2));
t1 = abs(pointlUp'*linerUpW)   ;
t2 = abs(pointlBot'*linerBotW) ;
t3 = abs(pointrUp'*linelUpW)   ;
t4 = abs(pointrBot'*linelBotW);
t5 = abs(pointUp'*ls);
t6 = abs(pointBot'*ls);
cost = [cost t1 t2 t3 t4 t5 t6];

% figure(6)
% hold off
% ShowPoly( data.Ov21, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
% hold on
% axis ij
% axis equal
% ShowPoly( data.Ov12, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15); 
% drawLine(linelUp,'r-');
% drawLine(linelBot,'r-');
% drawLine(linerUp,'b-');
% drawLine(linerBot,'b-');
% plot(ev1(1),ev1(2),'k+');
% plot(ev2(1),ev2(2),'k+');
% drawLine(ls,'y-');
% drawLine(linerUpW,'g--');
% drawLine(linerBotW,'g--');
% drawLine(linelUpW,'g--');
% drawLine(linelBotW,'g--');


return;
