function [cost ls] = getvy(a, b, polyBoundaryVecCell, data, vz, T)
data.vz = vz;
epsilon=1e-9;
[cost vy ls] = Golden(a, b, polyBoundaryVecCell, data, T, epsilon);

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
lh = data.lh;
fp = data.fp;
l=[-lh(2);lh(1);lh(2)*fp(1)-lh(1)*fp(2)];
perpLine = -sign(l(3))*l/norm(l(1:2));
fig=figure(5);
ShowPoly( Or, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
hold on
axis ij;
axis equal;
ShowPoly( Ov1, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15); 
ShowPoly( Ov2, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
ShowPoly( Ov12, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
ShowPoly( Ov21, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15);  
drawLine(lh, 'g-');
drawLine(perpLine, 'g-');
plot(data.fp(1),data.fp(2),'or');
plot(vy(1),vy(2),'r*');
close(fig);

function [cost vy ls] = Golden(a, b, polyBoundaryVecCell, data, T, epsilon)
lh = data.lh;
fp = data.fp;
l=[-lh(2);lh(1);lh(2)*fp(1)-lh(1)*fp(2)];
perpLine = -sign(l(3))*l/norm(l(1:2));

vy2(1,1)=a(1)+0.618*(b(1)-a(1));
vy2(2,1)=(-perpLine(1)*vy2(1)-perpLine(3))/perpLine(2);
[cost2 ~] = getvyls(polyBoundaryVecCell, data, vy2, T);
vy1(1,1)=a(1)+0.382*(b(1)-a(1));
vy1(2,1)=(-perpLine(1)*vy1(1)-perpLine(3))/perpLine(2);
[cost1 ~] = getvyls(polyBoundaryVecCell, data, vy1, T);
while(distPnts(a,b)>epsilon)
    if cost1<cost2
        b=vy2;
        vy2=vy1;
        cost2=cost1;
        vy1(1,1)=a(1)+0.382*(b(1)-a(1));
        vy1(2,1)=(-perpLine(1)*vy1(1)-perpLine(3))/perpLine(2);
        [cost1 ~] = getvyls(polyBoundaryVecCell, data, vy1, T);
        elseif cost1==cost2
            a=vy1;
            b=vy2;
            vy2(1,1)=a(1)+0.618*(b(1)-a(1));
            vy2(2,1)=(-perpLine(1)*vy2(1)-perpLine(3))/perpLine(2);            
            [cost2 ~] = getvyls(polyBoundaryVecCell, data, vy2, T);            
            vy1(1,1)=a(1)+0.382*(b(1)-a(1));
            vy1(2,1)=(-perpLine(1)*vy1(1)-perpLine(3))/perpLine(2);
            [cost1 ~] = getvyls(polyBoundaryVecCell, data, vy1, T);
    else
        a=vy1;
        vy1=vy2;
        cost1=cost2;
        vy2(1,1)=a(1)+0.618*(b(1)-a(1));
        vy2(2,1)=(-perpLine(1)*vy2(1)-perpLine(3))/perpLine(2);
        [cost2 ~] = getvyls(polyBoundaryVecCell, data, vy2, T);    
    end
end
vy=(a+b)/2;
[cost ls] = getvyls(polyBoundaryVecCell, data, vy, T);

return;


