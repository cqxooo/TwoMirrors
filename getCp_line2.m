function [cp c] = getCp_line2(polyBoundaryVecCell, lh, v, T)

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

v1 = v(:,1);
v2 = v(:,2);
v3 = v(:,3);
v4 = v(:,4);
p1 = ones(3,1);
p2 = ones(3,1);
[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Or',v1);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(Ov1',v1);
if flagl == 0 || flagr == 0
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
pointlUp = pointlUp/pointlUp(3);
pointrUp = pointrUp/pointrUp(3);
pointlBot = pointlBot/pointlBot(3);
pointrBot = pointrBot/pointrBot(3);
l1 = [linelUp linerUp linelBot linerBot];
p1(2) = (v1(2)*(pointlUp(2)+pointrUp(2))-2*pointlUp(2)*pointrUp(2))/(2*v1(2)-pointlUp(2)-pointrUp(2));
p1(1) = (v1(1)*(pointlUp(1)+pointrUp(1))-2*pointlUp(1)*pointrUp(1))/(2*v1(1)-pointlUp(1)-pointrUp(1));

p2(2) = (v1(2)*(pointlBot(2)+pointrBot(2))-2*pointlBot(2)*pointrBot(2))/(2*v1(2)-pointlBot(2)-pointrBot(2));
p2(1) = (v1(1)*(pointlBot(1)+pointrBot(1))-2*pointlBot(1)*pointrBot(1))/(2*v1(1)-pointlBot(1)-pointrBot(1));
l = cross(p1, p2);
m1 = -sign(l(3))*l/norm(l(1:2));


[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Or',v2);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(Ov2',v2);
if flagl == 0 || flagr == 0
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
pointlUp = pointlUp/pointlUp(3);
pointrUp = pointrUp/pointrUp(3);
pointlBot = pointlBot/pointlBot(3);
pointrBot = pointrBot/pointrBot(3);
l2 = [linelUp linerUp linelBot linerBot];
p1(2) = (v2(2)*(pointlUp(2)+pointrUp(2))-2*pointlUp(2)*pointrUp(2))/(2*v2(2)-pointlUp(2)-pointrUp(2));
p1(1) = (v2(1)*(pointlUp(1)+pointrUp(1))-2*pointlUp(1)*pointrUp(1))/(2*v2(1)-pointlUp(1)-pointrUp(1));

p2(2) = (v2(2)*(pointlBot(2)+pointrBot(2))-2*pointlBot(2)*pointrBot(2))/(2*v2(2)-pointlBot(2)-pointrBot(2));
p2(1) = (v2(1)*(pointlBot(1)+pointrBot(1))-2*pointlBot(1)*pointrBot(1))/(2*v2(1)-pointlBot(1)-pointrBot(1));
l = cross(p1, p2);
m2 = -sign(l(3))*l/norm(l(1:2));

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Ov1',v3);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(Ov12',v3);
if flagl == 0 || flagr == 0
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
pointlUp = pointlUp/pointlUp(3);
pointrUp = pointrUp/pointrUp(3);
pointlBot = pointlBot/pointlBot(3);
pointrBot = pointrBot/pointrBot(3);
l3 = [linelUp linerUp linelBot linerBot];
p1(2) = (v3(2)*(pointlUp(2)+pointrUp(2))-2*pointlUp(2)*pointrUp(2))/(2*v3(2)-pointlUp(2)-pointrUp(2));
p1(1) = (v3(1)*(pointlUp(1)+pointrUp(1))-2*pointlUp(1)*pointrUp(1))/(2*v3(1)-pointlUp(1)-pointrUp(1));

p2(2) = (v3(2)*(pointlBot(2)+pointrBot(2))-2*pointlBot(2)*pointrBot(2))/(2*v3(2)-pointlBot(2)-pointrBot(2));
p2(1) = (v3(1)*(pointlBot(1)+pointrBot(1))-2*pointlBot(1)*pointrBot(1))/(2*v3(1)-pointlBot(1)-pointrBot(1));
l = cross(p1, p2);
m3 = -sign(l(3))*l/norm(l(1:2));

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Ov21',v4);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(Ov2',v4);
if flagl == 0 || flagr == 0
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
pointlUp = pointlUp/pointlUp(3);
pointrUp = pointrUp/pointrUp(3);
pointlBot = pointlBot/pointlBot(3);
pointrBot = pointrBot/pointrBot(3);
l4 = [linelUp linerUp linelBot linerBot];
p1(2) = (v4(2)*(pointlUp(2)+pointrUp(2))-2*pointlUp(2)*pointrUp(2))/(2*v4(2)-pointlUp(2)-pointrUp(2));
p1(1) = (v4(1)*(pointlUp(1)+pointrUp(1))-2*pointlUp(1)*pointrUp(1))/(2*v4(1)-pointlUp(1)-pointrUp(1));

p2(2) = (v4(2)*(pointlBot(2)+pointrBot(2))-2*pointlBot(2)*pointrBot(2))/(2*v4(2)-pointlBot(2)-pointrBot(2));
p2(1) = (v4(1)*(pointlBot(1)+pointrBot(1))-2*pointlBot(1)*pointrBot(1))/(2*v4(1)-pointlBot(1)-pointrBot(1));
l = cross(p1, p2);
m4 = -sign(l(3))*l/norm(l(1:2));

% figure(6)
% hold off
% ShowPoly( Ov21, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
% hold on
% axis ij
% axis equal
% ShowPoly( Ov2, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15); 
% drawLine(linelUp,'r-');
% drawLine(linelBot,'r-');
% drawLine(linerUp,'b-');
% drawLine(linerBot,'b-');
% plot(pointlUp(1),pointlUp(2),'r+');
% plot(pointrUp(1),pointrUp(2),'r+');
% plot(pointlBot(1),pointlBot(2),'k+');
% plot(pointrBot(1),pointrBot(2),'k+');
% drawLine(mu4,'k-');
% drawLine(mb4,'k-');

lh = -sign(lh(3))*lh/norm(lh(1:2));
P = [1 0 0; 0 1 0; lh'];
l = inv(P)' * [l1 l2 l3 l4];
m = inv(P)' * [m1 m2 m3 m4];

Or_ = P*Or;
Or_ = [Or_(1,:)./Or_(3,:); Or_(2,:)./Or_(3,:); Or_(3,:)./Or_(3,:)];
Ov1_ = P*Ov1;
Ov1_ = [Ov1_(1,:)./Ov1_(3,:); Ov1_(2,:)./Ov1_(3,:); Ov1_(3,:)./Ov1_(3,:)];
Ov2_ = P*Ov2;
Ov2_ = [Ov2_(1,:)./Ov2_(3,:); Ov2_(2,:)./Ov2_(3,:); Ov2_(3,:)./Ov2_(3,:)];
Ov12_ = P*Ov12;
Ov12_ = [Ov12_(1,:)./Ov12_(3,:); Ov12_(2,:)./Ov12_(3,:); Ov12_(3,:)./Ov12_(3,:)];
Ov21_ = P*Ov21;
Ov21_ = [Ov21_(1,:)./Ov21_(3,:); Ov21_(2,:)./Ov21_(3,:); Ov21_(3,:)./Ov21_(3,:)];
fig=figure(4);
hold off;
ShowPoly( Or_, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
hold on
axis ij;   
axis equal;
ShowPoly( Ov1_, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15);  
ShowPoly( Ov2_, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
ShowPoly( Ov12_, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
ShowPoly( Ov21_, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15);   
for i=1:size(l,2)
    s = l(:,i);
    l(:,i) = -sign(s(3))*s/norm(s(1:2));
    drawLine(l(:,i),'g-');       
end

for i=1:size(m,2)
    s = m(:,i);
    m(:,i) = -sign(s(3))*s/norm(s(1:2));
    drawLine(m(:,i),'r-');  
end

l1 = l(1:2,1);
l2 = l(1:2,5);
l3 = l(1:2,9);
l4 = l(1:2,13);

m1 = m(1:2,1);
m2 = m(1:2,2);
m3 = m(1:2,3);
m4 = m(1:2,4);


lm1 = [l1 m1];
x=-lm1(2,:)./lm1(1,:);
lm2 = [m2 l2];
y=-lm2(2,:)./lm2(1,:);

a1 = x(1); b1 = x(2); a2 = y(1); b2 = y(2);

cen1 = (a1*b2-a2*b1)/(a1-b1-a2+b2);
r1 = sqrt(cen1^2+(a1-b1)*(a1*b1-a2*b2)/(a1-b1-a2+b2)-a1*b1);
c1=[1 0 -cen1; 0 1 0; -cen1 0 cen1^2-r1^2];
fig=figure(4);
drawConic(c1,'b-');
hold on;


lm3 = [l3 m3];
x=-lm3(2,:)./lm3(1,:);
lm4 = [m4 l4];
y=-lm4(2,:)./lm4(1,:);
a1 = x(1); b1 = x(2); a2 = y(1); b2 = y(2);
cen2 = (a1*b2-a2*b1)/(a1-b1-a2+b2);
r2 = sqrt(cen2^2+(a1-b1)*(a1*b1-a2*b2)/(a1-b1-a2+b2)-a1*b1);
c2=[1 0 -cen2; 0 1 0; -cen2 0 cen2^2-r2^2];
fig=figure(4);
drawConic(c2,'b-');

alpha = ((cen1^2-r1^2)-(cen2^2-r2^2))/(2*(cen1-cen2));
beta = sqrt(r1^2-(alpha-cen1)^2);

c = [(alpha-sqrt(-1)*beta)*lh(3); lh(3); -alpha*lh(1)-lh(2)+sqrt(-1)*beta*lh(1)];
c = [conj(c) c];
cp = [c(1,:)./c(3,:);c(2,:)./c(3,:);c(3,:)./c(3,:)];
close(fig);






