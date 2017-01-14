function [RR alpha beta] = rectify(polyBoundaryVecCell, K, ls, lh)
%===========================================
%rectify
%===========================================
x0 = K(1:2,3);
fl = K(1,1);
pbv=polyBoundaryVecCell{1};
s = size(pbv,2);
Or =  [pbv; ones(1,s)];
pbv=polyBoundaryVecCell{2};
s = size(pbv,2);
Ov1 = [pbv; ones(1,s)];
pbv=polyBoundaryVecCell{3};
s = size(pbv,2);
Ov2 = [pbv; ones(1,s)];
pbv=polyBoundaryVecCell{4};
s = size(pbv,2);
Ov12 =  [pbv; ones(1,s)];
pbv=polyBoundaryVecCell{5};
s = size(pbv,2);
Ov21 =  [pbv; ones(1,s)];

fig = 0;

% x_ = getPerPtLine(x0, ls);
% v1 = inv(K)*[x0;1]; v1=v1/norm(v1);
% v2 = inv(K)*[x_;1]; v2=v2/norm(v2);
% ct = v1'*v2;
% theta = acos(ct);
% v = cross(v2,v1);
% v = v/norm(v);
% R0 = rot(theta,v,0);
x_ = (-ls(2)*x0(2)-ls(3))/ls(1);
d = x0(1) - x_;
theta = atan2(d,fl);
R0 = rot(theta,'y',0);
H0 = K*R0*inv(K);

ls0 = inv(H0')*ls;
ls0 = -sign(ls0(3))*ls0/norm(ls0(1:2));

if fig
    fig1=figure(99);
    H = H0;    
    Or_ = H*Or;   
    Ov1_ = H*Ov1;   
    Ov2_ = H*Ov2;   
    Ov12_ = H*Ov12;   
    Ov21_ = H*Ov21;   
    Or_ = [Or_(1,:)./Or_(3,:);Or_(2,:)./Or_(3,:);Or_(3,:)./Or_(3,:)];
    Ov1_ = [Ov1_(1,:)./Ov1_(3,:);Ov1_(2,:)./Ov1_(3,:);Ov1_(3,:)./Ov1_(3,:)];
    Ov2_ = [Ov2_(1,:)./Ov2_(3,:);Ov2_(2,:)./Ov2_(3,:);Ov2_(3,:)./Ov2_(3,:)];
    Ov12_ = [Ov12_(1,:)./Ov12_(3,:);Ov12_(2,:)./Ov12_(3,:);Ov12_(3,:)./Ov12_(3,:)];
    Ov21_ = [Ov21_(1,:)./Ov21_(3,:);Ov21_(2,:)./Ov21_(3,:);Ov21_(3,:)./Ov21_(3,:)];
    ShowPoly( Or_, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
    hold on
    ShowPoly( Ov1_, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15); 
    ShowPoly( Ov2_, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
    ShowPoly( Ov12_, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
    ShowPoly( Ov21_, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15); 
    axis ij;
    plot(x0(1), x0(2),'r*','markersize',12);
    plot(x_, x0(2),'b*','markersize',12);
    drawLine(inv(H)'*ls,'b-');
    drawLine(inv(H)'*lh,'g-');    
end    
   
%     axis([0 640 0 480]);
%     set(gca, 'xtick', [], 'ytick', []);


%===========================================
% rotate the camera so that the axis is vertical
R1 = [ ls0(1) ls0(2) 0
     -ls0(2) ls0(1) 0
    0      0    1];

% theta = pi/2-atan(-ls0(1)/ls0(2));
% R1 = rot(theta,'z',0);

H1 = K*R1*R0*inv(K);
Or_ = H1*Or;   
Or_ = [Or_(1,:)./Or_(3,:); Or_(2,:)./Or_(3,:); Or_(3,:)./Or_(3,:)];
vz = cross(inv(H1)'*ls, inv(H1)'*lh);    
vz = vz./vz(3);
if Or_(2,:)<vz(2)    
    R1 = [ -ls0(1) -ls0(2) 0
            ls0(2) -ls0(1) 0
               0      0    1];
    H1 = K*R1*R0*inv(K);
end


%===========================================
% rectified contur in image

if fig
    fig2=figure(98);
    H = H1;
    Or_ = H*Or;   
    Ov1_ = H*Ov1;   
    Ov2_ = H*Ov2;   
    Ov12_ = H*Ov12;   
    Ov21_ = H*Ov21;   
    Or_ = [Or_(1,:)./Or_(3,:);Or_(2,:)./Or_(3,:);Or_(3,:)./Or_(3,:)];
    Ov1_ = [Ov1_(1,:)./Ov1_(3,:);Ov1_(2,:)./Ov1_(3,:);Ov1_(3,:)./Ov1_(3,:)];
    Ov2_ = [Ov2_(1,:)./Ov2_(3,:);Ov2_(2,:)./Ov2_(3,:);Ov2_(3,:)./Ov2_(3,:)];
    Ov12_ = [Ov12_(1,:)./Ov12_(3,:);Ov12_(2,:)./Ov12_(3,:);Ov12_(3,:)./Ov12_(3,:)];
    Ov21_ = [Ov21_(1,:)./Ov21_(3,:);Ov21_(2,:)./Ov21_(3,:);Ov21_(3,:)./Ov21_(3,:)];
    ShowPoly( Or_, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
    hold on
    ShowPoly( Ov1_, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15); 
    ShowPoly( Ov2_, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
    ShowPoly( Ov12_, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
    ShowPoly( Ov21_, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15); 
    axis ij; 
    drawLine(inv(H)'*ls,'b-');
    drawLine(inv(H)'*lh,'g-');       
end

%     axis([0 640 0 480]);
%     set(gca, 'xtick', [], 'ytick', []);



%===========================================
% rotate the camera so that the intersection of ls and lh align with
% principle point.

% x_ = getPerPtLine(x0, inv(H1)'*lh);
% v1 = inv(K)*[x0;1]; v1=v1/norm(v1);
% v2 = inv(K)*[x_;1]; v2=v2/norm(v2);
% ct = v1'*v2;
% theta = acos(ct);
% v = cross(v2,v1);
% v = v/norm(v);
% R2 = rot(theta,v,0);
lh0 = inv(H1')*lh;
lh0 = -sign(lh0(3))*lh0/norm(lh0(1:2));
y_ = (-lh0(1)*x0(1)-lh0(3))/lh0(2);
d = y_ - x0(2);
theta = atan2(d,fl);
R2 = rot(theta,'x',0);
H2 = K*R2*R1*R0*inv(K);

ls0 = inv(H2')*ls;
ls0 = -sign(ls0(3))*ls0/norm(ls0(1:2))

lh0 = inv(H2')*lh;
lh0 = -sign(lh0(3))*lh0/norm(lh0(1:2))

if fig
    fig3=figure(97);
    H = H2;
    Or_ = H*Or;   
    Ov1_ = H*Ov1;   
    Ov2_ = H*Ov2;   
    Ov12_ = H*Ov12;   
    Ov21_ = H*Ov21;   
    Or_ = [Or_(1,:)./Or_(3,:);Or_(2,:)./Or_(3,:);Or_(3,:)./Or_(3,:)];
    Ov1_ = [Ov1_(1,:)./Ov1_(3,:);Ov1_(2,:)./Ov1_(3,:);Ov1_(3,:)./Ov1_(3,:)];
    Ov2_ = [Ov2_(1,:)./Ov2_(3,:);Ov2_(2,:)./Ov2_(3,:);Ov2_(3,:)./Ov2_(3,:)];
    Ov12_ = [Ov12_(1,:)./Ov12_(3,:);Ov12_(2,:)./Ov12_(3,:);Ov12_(3,:)./Ov12_(3,:)];
    Ov21_ = [Ov21_(1,:)./Ov21_(3,:);Ov21_(2,:)./Ov21_(3,:);Ov21_(3,:)./Ov21_(3,:)];
    ShowPoly( Or_, 'FaceColor', MyPalette(1),  'EdgeColor',  MyPalette(1), 'FaceAlpha', 0.15);  
    hold on
    ShowPoly( Ov1_, 'FaceColor', MyPalette(2),  'EdgeColor',  MyPalette(2), 'FaceAlpha', 0.15); 
    ShowPoly( Ov2_, 'FaceColor', MyPalette(3),  'EdgeColor',  MyPalette(3), 'FaceAlpha', 0.15);  
    ShowPoly( Ov12_, 'FaceColor', MyPalette(4),  'EdgeColor',  MyPalette(4), 'FaceAlpha', 0.15);  
    ShowPoly( Ov21_, 'FaceColor', MyPalette(5),  'EdgeColor',  MyPalette(5), 'FaceAlpha', 0.15); 
    axis ij;  
    plot(x0(1), x0(2),'r*','markersize',12);
    plot(x0(1), y_,'b*','markersize',12);
    drawLine(inv(H)'*ls,'b-');
    drawLine(inv(H)'*lh,'g-'); 
end

RR = inv(R2*R1*R0);
if fig
    close(fig1);
    close(fig2);
    close(fig3);
end
% v = H2*vp;
% v = [v(1,:)./v(3,:); v(2,:)./v(3,:); v(3,:)./v(3,:)];
% plot(v(1,1),v(2,1),'r*','MarkerSize',12);
% plot(v(1,2),v(2,2),'b*','MarkerSize',12);
% cp = H2*cp;
% cp = [cp(1,:)./cp(3,:); cp(2,:)./cp(3,:); cp(3,:)./cp(3,:)];
% [alpha beta] = getAngle(v,inv(H2)'*lh, inv(H2)'*ls, cp);

return;