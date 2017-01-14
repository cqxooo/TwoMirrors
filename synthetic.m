close all
clear all
clc


lss = [];
lhs = [];
cps = [];
vxs = [];
thetas = [];
data={};
K9 = {};
angle9 = [];
R = {};
PrjMat = {};
cen = {};
x = {};
conic = {};
polyBoundaryVecCell = {};
imgSize = [640,480];

f       = 1500;
s = 0;
alpha   = 1;
pp = [imgSize(1), imgSize(2)]/2;
K0 = [alpha*f s pp(1);
    0        f pp(2);
    0        0  1    ];

%===1 trial
% R =  [ 0.8648    0.3161   -0.3901;
%       -0.4936    0.6778   -0.5450;
%        0.0921    0.6639    0.7421];

%===2 trial
% R{2} = rot(30, 'x', 1)*rot(-30, 'z', 1)*rot(-30, 'y', 1);
% %===3 trial
% R{3} =  [ 0.8648      0    -0.3901;
%           0        1        0  ;
%        0.0921      0     0.7421]*rot(30, 'x', 1)*rot(-30, 'z', 1)*rot(-30, 'y', 1);
% %===4 trial
% R{4} =  [    1       0         0  ;
%           0    0.6778   -0.5450;
%           0    0.6639    0.7421]*rot(30, 'x', 1)*rot(-30, 'z', 1)*rot(-30, 'y', 1);
% % ===5 trial
% R{5} =  [ 0.8648    0.3161       0;
%       -0.4936    0.6778       0;
%         0          0          1]*rot(30, 'x', 1);
t = [0 -170 -150]'; 
sigma = [-1 0 0;0 1 0;0 0 1];
alpha = 30*pi/180;
beta =  42*pi/180;

r = 10;  
cen{1}=[0 -10 -30]'; 
cen{2} = rot(2*alpha, 'y', 0)*sigma*cen{1};
cen{3} = rot(2*pi-2*beta, 'y', 0)*sigma*cen{1};
cen{4} = rot(2*alpha+2*beta, 'y', 0)*cen{1};
cen{5} = rot(2*pi-2*alpha-2*beta, 'y', 0)*cen{1};
fig = figure(1);
% R = rot(10, 'x', 1)*rot(-30, 'z', 1)*rot(10, 'y', 1);
R = rot(0, 'x', 1)*rot(-5, 'z', 1)*rot(5, 'y', 1);

PrjMat = K0*R*[eye(3) -t];
cp0 = [];
cp = PrjMat*[1 0 sqrt(-1) 0]';%cp circule points
cp0 = [cp0 cp./cp(3)]; 
cp = PrjMat*[1 0 -sqrt(-1) 0]';
cp0 = [cp0 cp./cp(3)]; 
v = PrjMat*[0 0 1 0]';
vz0 = v./v(3);
v = PrjMat*[0 1 0 0]';
vy0 = v./v(3);
v = PrjMat*[1 0 0 0]';
vx0 = v./v(3);
l = cross(cp0(:,1),cp0(:,2));
lh0 = -sign(l(3))*l/norm(l(1:2));
l = cross(vz0,vy0);
ls0 = -sign(l(3))*l/norm(l(1:2));
fig = figure(1);
    for j = 1:5        
        Q{j} = [eye(3) -cen{j}; -cen{j}' cen{j}(1)^2+cen{j}(2)^2+cen{j}(3)^2-r^2];
        C = inv(PrjMat*inv(Q{j})*PrjMat');
        hold on
        [~, x{j}] = drawConic(C,'b-');
        polyBoundaryVec{j} = x{j}(:,1:2)';
        conic{j} = C;
        axis ij, axis equal, axis equal
    end
    drawLine(lh0, 'b-');
    drawLine(ls0, 'b-');
 close(fig);


% c1 = inv(C1(1:2,1:2))*(-C1(1:2,3));
% c2 = inv(C2(1:2,1:2))*(-C2(1:2,3));
% p = intersectConics(C1, C2);
for sigma = 0:10/20:10
    sigma 
    num = 100;
for imgLoop = 1:num    
    imgLoop
    polyBoundaryVecCell{imgLoop} = polyBoundaryVec;    

    if sigma~=0
    for viewLoop = 1:length(polyBoundaryVecCell{imgLoop})
        N = size(polyBoundaryVecCell{imgLoop}{viewLoop},2);
%         while 1
            e = -1*sigma + 2*sigma * randn(1,N);
%             if max(e)<sigma && min(e)>-sigma
%                 break;
%             end
%         end
        for i = 1:N
            cur_1=mod(i-1+N-1,N)+1;
            cur_2=mod(i-2+N-1,N)+1;
            cur1=mod(i+1+N-1,N)+1;
            cur2=mod(i+2+N-1,N)+1;
            t = [e(i) e(cur_1) e(cur1) e(cur_2) e(cur2)];
            %e(i) = median(t);
            te=4*e(i)+2*e(cur_1)+2*e(cur1)+e(cur_2)+e(cur2);
            e(i)=te/10.0;
        end
        e = [e;e;ones(1,N)];
        for i=1:size(polyBoundaryVecCell{imgLoop}{viewLoop},2)
            s = [polyBoundaryVecCell{imgLoop}{viewLoop};ones(1,size(polyBoundaryVecCell{imgLoop}{viewLoop},2))];
            l = conic{viewLoop}*s(:,i);
            l = -sign(l(3))*l/norm(l(1:2));
            d(:,i)=[l(1:2); 0];
        end
        s = s+e.*d;
        s = [s(1,:)./s(3,:); s(2,:)./s(3,:); s(3,:)./s(3,:)];
        polyBoundaryVecCell{imgLoop}{viewLoop} = s(1:2,:);
    end
    end
    fig = figure(1);     
        for viewLoop = 1:5,
            hold on
           ShowPoly( polyBoundaryVecCell{imgLoop}{viewLoop},...
                'FaceColor', MyPalette(viewLoop),  'EdgeColor',  MyPalette(viewLoop),      'FaceAlpha', 0.15);    
            axis ij, axis off            
        end
     close(fig);
     T=normalize(polyBoundaryVecCell{imgLoop});   
    
    [lh,v]=getLh(polyBoundaryVecCell{imgLoop}, T);   
    
    [lh1,v0,v1, residual]=getRefinedLh(polyBoundaryVecCell{imgLoop},lh,v,T);
    
    [cp c]  = getCp_angle(polyBoundaryVecCell{imgLoop},lh1, v1, T);
    
    theta = getTheta(cp,v1);
    
    EvalPrint('(theta) * (180/pi)');

    ls = Initvzls(polyBoundaryVecCell{imgLoop},cp,theta,lh1,v1,T);
    
%     ls = InitLs_new(polyBoundaryVecCell{imgLoop},cp,theta,lh1,v1,T);
    data{imgLoop}= getRefinedLs_new(polyBoundaryVecCell{imgLoop}, lh1, ls, cp, theta, T);

    
    thetas = [thetas theta];
    lss = [lss data{imgLoop}.ls];
    lhs = [lhs data{imgLoop}.lh];
    cps = [cps data{imgLoop}.cp];
    vxs = [vxs data{imgLoop}.vx];
end
   [error, IAC, K] = calib_cplsvxs(cps, lss, vxs, K0);
   eval(['K' int2str(2*sigma+1) '=K;']);
   eval(['thetas' int2str(2*sigma+1) '=mean(thetas);']);
   save(['syn' int2str(2*sigma+1) '.mat'],['K' int2str(2*sigma+1)],['thetas' int2str(2*sigma+1)]);
end
   