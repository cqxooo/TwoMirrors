function [cp c] = getCircle(polyBoundaryVecCell,lh, v, T)

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

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Or',v1);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(Ov1',v1);
if flagl == 0 || flagr == 0
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
l1 = [linelUp linelBot linerUp linerBot];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Or',v2);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(Ov2',v2);
if flagl == 0 || flagr == 0
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
l2 = [linelUp linelBot linerUp linerBot];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Ov1',v3);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(Ov12',v3);
if flagl == 0 || flagr == 0
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
l3 = [linelUp linelBot linerUp linerBot];

[linelUp,linelBot,pointlUp,pointlBot,flagl] = getOuterTengents(Ov2',v4);
[linerUp,linerBot,pointrUp,pointrBot,flagr] = getOuterTengents(Ov21',v4);
if flagl == 0 || flagr == 0
    return;
end
linelUp = -sign(linelUp(3))*linelUp/norm(linelUp(1:2));
linelBot = -sign(linelBot(3))*linelBot/norm(linelBot(1:2));
linerUp = -sign(linerUp(3))*linerUp/norm(linerUp(1:2));
linerBot = -sign(linerBot(3))*linerBot/norm(linerBot(1:2));
l4 = [linelUp linelBot linerUp linerBot];

% figure(6)
% hold off
% ShowPoly( Ov2, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
% hold on
% axis ij
% axis equal
% ShowPoly( Ov21, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15); 
% drawLine(linelUp,'r-');
% drawLine(linelBot,'r-');
% drawLine(linerUp,'b-');
% drawLine(linerBot,'b-');
% plot(v2(1),v2(2),'k+');
% drawLine(lh,'y-');

lh = -sign(lh(3))*lh/norm(lh(1:2));
P = [1 0 0; 0 1 0; lh'];
l = inv(P)'*[l1 l2 l3 l4];

fig=figure(4);
hold off;
ShowPoly( Or, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
hold on
axis ij;
axis equal;   
ShowPoly( Ov1, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15);  
ShowPoly( Ov2, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
ShowPoly( Ov12, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
ShowPoly( Ov21, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15);  

for i = 1:size(l,2);
        s = l(:,i);
        l(:,i) = -sign(s(3))*s/norm(s(1:2));
        if fig
            S = 'g';
            drawLine(l(:,i),S);     
            hold on;
        end
end
close(fig)
l1 = l(1:2,1);
l2 = l(1:2,5);
l3 = l(1:2,9);
l4 = l(1:2,13);

l12 = [l4 l2];
x=-l12(2,:)./l12(1,:);
l12 = [l1 l3];
y=-l12(2,:)./l12(1,:);

a1 = x(1); b1 = x(2); a2 = y(1); b2 = y(2);

cen1 = (a1*b2-a2*b1)/(a1-b1-a2+b2);
r1 = sqrt(cen1^2+(a1-b1)*(a1*b1-a2*b2)/(a1-b1-a2+b2)-a1*b1);
c1=[1 0 -cen1; 0 1 0; -cen1 0 cen1^2-r1^2];
fig=figure(4);
drawConic(c1,'b-');
hold on;


l12 = [l2 l1];
y=-l12(2,:)./l12(1,:);
a2 = y(1); b2 = y(2);
cen2 = (a1*b2-a2*b1)/(a1-b1-a2+b2);
r2 = sqrt(cen2^2+(a1-b1)*(a1*b1-a2*b2)/(a1-b1-a2+b2)-a1*b1);
c2=[1 0 -cen2; 0 1 0; -cen2 0 cen2^2-r2^2];
drawConic(c2,'b-');
close(fig);

alpha = ((cen1^2-r1^2)-(cen2^2-r2^2))/(2*(cen1-cen2));
beta = sqrt(r1^2-(alpha-cen1)^2);

c = [(alpha-sqrt(-1)*beta)*lh(3); lh(3); -alpha*lh(1)-lh(2)+sqrt(-1)*beta*lh(1)];
c = [conj(c) c];
cp = [c(1,:)./c(3,:);c(2,:)./c(3,:);c(3,:)./c(3,:)];

[lineUp,lineBot,pointlUp,pointrUp,pointlBot,pointrBot] = getBitangents(Or, Ov1);
% fig=figure(4);
% hold off;
% ShowPoly( Or, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
% hold on
% axis ij;
% axis equal;   
% ShowPoly( Ov1, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15);  
% plot(pointlUp(1),pointlUp(2), 'r+');
% plot(pointrUp(1),pointrUp(2), 'r+');
A = [1/beta -alpha/beta 0;0 1 0;0 0 1];
pl = A*P*pointlUp;
pl = pl/pl(3);
pr = A*P*pointrUp;
pr = pr/pr(3);
cc = A*P*c; 

syms a b r;
eqn1 = (pl(1)-a)^2+(pl(2)-b)^2-r^2;
eqn2 = (pr(1)-a)^2+(pr(2)-b)^2-r^2;
eqn3 = (cc(1,1)-a)^2+(cc(2,1)-b)^2-r^2;
eqn4 = (cc(1,2)-a)^2+(cc(2,2)-b)^2-r^2;
answer = solve(eqn1, eqn2, eqn3, eqn4, 'a','b','r');

% drawLine(ls,'g-');
% r=crossRatio(v1,v2,v3,v4);
% angle=acos(1/(2*sqrt(r)));


return;