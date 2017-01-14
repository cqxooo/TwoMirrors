function [lh,v]=getLh2(polyBoundaryVecCell,T) 
% this function tries to refine the coeficent ktan(theta/2).
% data structure has ls, vx,lh,Cl, Cr.
% Bspline = getMirrorSilhouettes_n(filename);
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

figure(3)
hold off
ShowPoly( Or, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
hold on
axis ij;
axis equal;   
ShowPoly( Ov1, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15);  
ShowPoly( Ov2, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
ShowPoly( Ov12, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
ShowPoly( Ov21, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15);  

length0 = size(Or,2);
length1 = size(Ov1,2);
length2 = size(Ov2,2);
length12 = size(Ov12,2);
length21 = size(Ov21,2);

l1 = [];
l2 = [];
l3 = [];
l4 = [];

a = [Or'; Ov1'];
N = size(a,1);
length = length0;
idx = convhull(a(:,1),a(:,2));
% plot(x(idx),y(idx),'r-',x,y,'b+')
num= size(idx,1);
for i = 1:num
    if idx(i) <= length && idx(cycle(i+1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i+1,num)),:))';
        drawLine(l,'r-');
        l1 = [l1 l];
    end
    if idx(i) <= length && idx(cycle(i-1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i-1,num)),:))';
        drawLine(l,'r-');
        l1 = [l1 l];
    end

end
a = [Ov2'; Ov12'];
N = size(a,1);
length = length2;
idx = convhull(a(:,1),a(:,2));
% plot(x(idx),y(idx),'r-',x,y,'b+')
num= size(idx,1);
for i = 1:num
    if idx(i) <= length && idx(cycle(i+1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i+1,num)),:))';
        drawLine(l,'r-');
        l1 = [l1 l];
    end
    if idx(i) <= length && idx(cycle(i-1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i-1,num)),:))';
        drawLine(l,'r-');
        l1 = [l1 l];
    end

end

p1 = intersection(l1);


a = [Or'; Ov2'];
N = size(a,1);
length = length0;
idx = convhull(a(:,1),a(:,2));
% plot(x(idx),y(idx),'r-',x,y,'b+')
num= size(idx,1);
for i = 1:num
    if idx(i) <= length && idx(cycle(i+1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i+1,num)),:))';
        drawLine(l,'r-');
        l2 = [l2 l];
    end
    if idx(i) <= length && idx(cycle(i-1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i-1,num)),:))';
        drawLine(l,'r-');
        l2 = [l2 l];
    end

end
a = [Ov1'; Ov21'];
N = size(a,1);
length = length1;
idx = convhull(a(:,1),a(:,2));
% plot(x(idx),y(idx),'r-',x,y,'b+')
num= size(idx,1);
for i = 1:num
    if idx(i) <= length && idx(cycle(i+1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i+1,num)),:))';
        drawLine(l,'r-');
        l2 = [l2 l];
    end
    if idx(i) <= length && idx(cycle(i-1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i-1,num)),:))';
        drawLine(l,'r-');
        l2 = [l2 l];
    end

end

p2 = intersection(l2);



a = [Ov1'; Ov12'];
N = size(a,1);
length = length1;
idx = convhull(a(:,1),a(:,2));
% plot(x(idx),y(idx),'r-',x,y,'b+')
num= size(idx,1);
for i = 1:num
    if idx(i) <= length && idx(cycle(i+1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i+1,num)),:))';
        drawLine(l,'r-');
        l3 = [l3 l];
    end
    if idx(i) <= length && idx(cycle(i-1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i-1,num)),:))';
        drawLine(l,'r-');
        l3 = [l3 l];
    end

end

p3 = intersection(l3);


a = [Ov2'; Ov21'];
N = size(a,1);
length = length2;
idx = convhull(a(:,1),a(:,2));
% plot(x(idx),y(idx),'r-',x,y,'b+')
num= size(idx,1);
for i = 1:num
    if idx(i) <= length && idx(cycle(i+1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i+1,num)),:))';
        drawLine(l,'r-');
        l4 = [l4 l];
    end
    if idx(i) <= length && idx(cycle(i-1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i-1,num)),:))';
        drawLine(l,'r-');
        l4 = [l4 l];
    end

end

p4 = intersection(l4);

p = [p1 p2 p3 p4];
[RANSAC_mask, l] = RANSAC(p(1:2,:), @TLS, 2, @Dist, 50, 2);
lh = -sign(l(3))*l/norm(l(1:2));

fig=figure(3);
hold off
ShowPoly( Or, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
hold on
axis ij;
axis equal;   
ShowPoly( Ov1, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15);  
ShowPoly( Ov2, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
ShowPoly( Ov12, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
ShowPoly( Ov21, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15); 
drawLine(lh,'g-');

close(fig);
% v = inv(data.T)*[v1 v2 v3 v4];
% v = [ v(1,:)./v(3,:); v(2,:)./v(3,:); v(3,:)./v(3,:)];
% l = data.T'*lh;
% lh = -sign(l(3))*l/norm(l(1:2));
% 
return;


