function ls = InitLs_new(polyBoundaryVecCell, cp, theta, lh, v, T)
v = [ v(1,:)./v(3,:); v(2,:)./v(3,:); v(3,:)./v(3,:)];
v1=v(2,1);
v2=v(2,2);
i=cp(2,1);
j=cp(2,2);
% b = v(:,3);
% a = v(:,4);
q1(2) = (2*i*j-v1*i-v1*j)/(i+j-2*v1);
q1(1) = -(lh(3)+lh(2)*q1(2))/lh(1);
q1(3) = 1;
q1 = q1';

q2(2) = (2*i*j-v2*i-v2*j)/(i+j-2*v2);
q2(1) = -(lh(3)+lh(2)*q2(2))/lh(1);
q2(3) = 1;
q2 = q2';
b = q1;
a = q2;
% mid = a+b/2;
epsilon=1e-9;
[cost vz ls] = Golden(a, b, polyBoundaryVecCell, cp, theta, lh, T, epsilon);

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

function [cost vz ls] = Golden(a, b, polyBoundaryVecCell, cp, theta, lh, T, epsilon)
vz2=a+0.618*(b-a);
[cost2 ~] = getLs_new(polyBoundaryVecCell, cp, theta, lh, vz2, T);
vz1=a+0.382*(b-a);
[cost1 ~] = getLs_new(polyBoundaryVecCell, cp, theta, lh, vz1, T);
while(distPnts(a,b)>epsilon)
    if cost1<cost2
        b=vz2;
        vz2=vz1;
        cost2=cost1;
        vz1=a+0.382*(b-a);
        [cost1 ~] = getLs_new(polyBoundaryVecCell, cp, theta, lh, vz1, T);
        elseif cost1==cost2
            a=vz1;
            b=vz2;
            vz2=a+0.618*(b-a);
            [cost2 ~] = getLs_new(polyBoundaryVecCell, cp, theta, lh, vz2, T);
            vz1=a+0.382*(b-a);
            [cost1 ~] = getLs_new(polyBoundaryVecCell, cp, theta, lh, vz1, T);
    else
        a=vz1;
        vz1=vz2;
        cost1=cost2;
        vz2=a+0.618*(b-a);
        [cost2 ~] = getLs_new(polyBoundaryVecCell, cp, theta, lh, vz2, T);        
    end
end
vz=(a+b)/2;
[cost ls] = getLs_new(polyBoundaryVecCell, cp, theta, lh, vz, T);

return;


