function [cp c] = getCp_line1(polyBoundaryVecCell,imF, lh, T)

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

img=imread(imF);
fig=figure(4);
image(img);
hold on
points=zeros(3,5);
temp = ones(3,1);
for i=1:5
    disp('Please click a point of five views...');
    [temp(1), temp(2)] = ginput(1);
    temp = T*temp;
    switch i
        case 1
            d = [];
            for j=1:size(Or,2)
                d = [d distPnts(temp,Or(:,j))];                
            end
            [mind ind] = min(d);
            point(:,i) = Or(:,ind);
            out = inv(T)*point(:,i);
            plot(out(1), out(2), 'r*');  
       case 2
            d = [];
            for j=1:size(Ov1,2)
                d = [d distPnts(temp,Ov1(:,j))];                
            end
            [mind ind] = min(d);
            point(:,i) = Ov1(:,ind);
            out = inv(T)*point(:,i);
            plot(out(1), out(2), 'r*');  
       case 3
            d = [];
            for j=1:size(Ov2,2)
                d = [d distPnts(temp,Ov2(:,j))];                
            end
            [mind ind] = min(d);
            point(:,i) = Ov2(:,ind);
            out = inv(T)*point(:,i);
            plot(out(1), out(2), 'r*');    
       case 4
            d = [];
            for j=1:size(Ov12,2)
                d = [d distPnts(temp,Ov12(:,j))];                
            end
            [mind ind] = min(d);
            point(:,i) = Ov12(:,ind);
            out = inv(T)*point(:,i);
            plot(out(1), out(2), 'r*');    
       case 5
            d = [];
            for j=1:size(Ov21,2)
                d = [d distPnts(temp,Ov21(:,j))];                
            end
            [mind ind] = min(d);
            point(:,i) = Ov21(:,ind);
            out = inv(T)*point(:,i);
            plot(out(1), out(2), 'r*');             
    end    
end
close(fig);

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
if fig 
    for i = 1:size(point,2)
        plot(point(1,i),point(2,i),'r+');
        hold on;
    end
end
close(fig)
lh = -sign(lh(3))*lh/norm(lh(1:2));
P = [1 0 0; 0 1 0; lh'];
pc_ = P * point;
pc_ = [pc_(1,:)./pc_(3,:); pc_(2,:)./pc_(3,:); pc_(3,:)./pc_(3,:)];
l = cross(pc_(:,1), pc_(:,2));
l12 = -sign(l(3))*l/norm(l(1:2));
mp12 = (pc_(:,1)+pc_(:,2))/2;
l = [l12(2); -l12(1); -l12(2)*mp12(1)+l12(1)*mp12(2)];
pl12 = -sign(l(3))*l/norm(l(1:2));
l = cross(pc_(:,1), pc_(:,3));
l13 = -sign(l(3))*l/norm(l(1:2));
mp13 = (pc_(:,1)+pc_(:,3))/2;
l = [l13(2); -l13(1); -l13(2)*mp13(1)+l13(1)*mp13(2)];
pl13 = -sign(l(3))*l/norm(l(1:2));
l = cross(pc_(:,2), pc_(:,4));
l24 = -sign(l(3))*l/norm(l(1:2));
mp24 = (pc_(:,2)+pc_(:,4))/2;
l = [l24(2); -l24(1); -l24(2)*mp24(1)+l24(1)*mp24(2)];
pl24 = -sign(l(3))*l/norm(l(1:2));
l = cross(pc_(:,3), pc_(:,5));
l35 = -sign(l(3))*l/norm(l(1:2));
mp35 = (pc_(:,3)+pc_(:,5))/2;
l = [l35(2); -l35(1); -l35(2)*mp35(1)+l35(1)*mp35(2)];
pl35 = -sign(l(3))*l/norm(l(1:2));
l = cross(pc_(:,4), pc_(:,5));
l45 = -sign(l(3))*l/norm(l(1:2));
mp45 = (pc_(:,4)+pc_(:,5))/2;
l = [l45(2); -l45(1); -l45(2)*mp45(1)+l45(1)*mp45(2)];
pl45 = -sign(l(3))*l/norm(l(1:2));
A = [pl12 pl13 pl24 pl45 pl35];
[u d v] = svd(A');
cen = v(:,3);
cen = cen/cen(3);

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
axis equal;   
ShowPoly( Ov1_, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15);  
ShowPoly( Ov2_, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
ShowPoly( Ov12_, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
ShowPoly( Ov21_, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15);  
if fig 
    for i = 1:size(pc_,2)
        plot(pc_(1,i),pc_(2,i),'r+');
        hold on;
    end
end
drawLine(pl12, 'g-');
drawLine(pl13, 'g-');
drawLine(pl24, 'g-');
drawLine(pl45, 'g-');
drawLine(pl35, 'g-');
plot(cen(1),cen(2),'k*');

p1 = pc_(:,4);
p2 = pc_(:,2);
x1 = p1(1)-cen(1);
y1 = p1(2)-cen(2);
x2 = p2(1)-cen(1);
y2 = p2(2)-cen(2);
s = 1;
cen1 = (x1*y1-s^2*x2*y2)/(y1^2-s^2*y2^2);
r1 = abs(s*(x2*y1-x1*y2)/(y1^2-s^2*y2^2));
c1=[1 0 -cen1; 0 1 0; -cen1 0 cen1^2-r1^2];
fig=figure(4);
drawConic(c1,'b-');
hold on;


p3 = pc_(:,5);
x3 = p3(1)-cen(1);
y3 = p3(2)-cen(2);
s = 1;
cen2 = (x1*y1-s^2*x3*y3)/(y1^2-s^2*y3^2);
r2 = abs(s*(x3*y1-x1*y3)/(y1^2-s^2*y3^2));
c2=[1 0 -cen2; 0 1 0; -cen2 0 cen2^2-r2^2];
fig=figure(4);
drawConic(c2,'b-');

alpha = ((cen1^2-r1^2)-(cen2^2-r2^2))/(2*(cen1-cen2));
beta = sqrt(r1^2-(alpha-cen1)^2);

c = [(alpha-sqrt(-1)*beta)*lh(3); lh(3); -alpha*lh(1)-lh(2)+sqrt(-1)*beta*lh(1)];
c = [conj(c) c];
cp = [c(1,:)./c(3,:);c(2,:)./c(3,:);c(3,:)./c(3,:)];
close(fig);




return;