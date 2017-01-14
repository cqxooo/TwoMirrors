close all
clear all
clc

imgFilename = [];
lss = [];
lhs = [];
cps = [];
vxs = [];
thetas = [];
kappas = [];
angles = {};
v1s = {};
data={};

% n = 1;
% for i=1:4
%     eval([' imgFilename{end+1} =''robotics/robot' int2str(n) '/frame' int2str(i) '.jpg'';']);
% end
% n = 2;
% for i=1:4
%     eval([' imgFilename{end+1} =''robotics/robot' int2str(n) '/frame' int2str(i) '.jpg'';']);
% end
% n = 3;
% for i=1:4
%     eval([' imgFilename{end+1} =''robotics/robot' int2str(n) '/frame' int2str(i) '.jpg'';']);
% end
% for i=1:3
%     eval(['imgFilename{end+1} =''model/baymax/baymax' int2str(i)  '.jpg'';']);
% end
% for i=1:3
%     eval(['imgFilename{end+1} =''model/elephant/elephant' int2str(i)  '.jpg'';']);
% end
% for i=1:3
%     eval(['imgFilename{end+1} =''model/vase/vase' int2str(i)  '.jpg'';']);
% end
for i=1:3
    eval(['imgFilename{end+1} =''formike/girl/girl' int2str(i)  '.jpg'';']);
end


[pathstr,imNameC,imext] = fileparts(imgFilename{1}); 

% for imgLoop = 1:length(imgFilename)     
%     disp(imgFilename{imgLoop})  
%     imCell{imgLoop} = imread(imgFilename{imgLoop});
%     polyBoundaryVecCell{imgLoop} = loadFromTxt(imgFilename{imgLoop},5);
%     
%      figure
%         for viewLoop = 1:5,       
%             hold on
%             ShowPoly( polyBoundaryVecCell{imgLoop}{viewLoop},...
%                 'FaceColor', MyPalette(viewLoop),  'EdgeColor',  MyPalette(viewLoop),      'FaceAlpha', 0.15);    
%             axis ij,    axis equal               
%         end
%             
%     T=normalize(polyBoundaryVecCell{imgLoop});   
%     
%     [lh,v]=getLh(polyBoundaryVecCell{imgLoop}, T);   
%     
%     [lh1,v0,v1, residual]=getRefinedLh(polyBoundaryVecCell{imgLoop},lh,v,T);
% 
%     [cp c]  = getCp_angle(polyBoundaryVecCell{imgLoop},lh1, v1, T);
% %     [cp c]  = getCircle(polyBoundaryVecCell{imgLoop},lh1, v1, T);
% 
%     theta = getTheta(cp,v1);
%     
%     EvalPrint('(theta) * (180/pi)');
% % save test.mat
% % load test.mat
% 
%     ls = InitLs_new(polyBoundaryVecCell{imgLoop},cp,theta,lh1,v1,T);
% %     [vz ls] = Initvzls(polyBoundaryVecCell{imgLoop},cp,theta,lh1,v1,T);
% 
%     data{imgLoop}= getRefinedLs_new(polyBoundaryVecCell{imgLoop}, lh1, ls, cp, theta, T);
%     
%     figure(99)
%     img=imread(imgFilename{imgLoop});
%     image(img);
%     hold on
%     drawLine(data{imgLoop}.ls, 'g-');
%     drawLine(data{imgLoop}.lh, 'g-');
%     
%     [alpha beta]=getAngle(v1,data{imgLoop}.lh,data{imgLoop}.ls,data{imgLoop}.cp, T);
%     lss = [lss data{imgLoop}.ls];
%     lhs = [lhs data{imgLoop}.lh];
%     cps = [cps data{imgLoop}.cp];
%     vxs = [vxs data{imgLoop}.vx];
%     angles{imgLoop} = [alpha beta];
%     thetas = [thetas theta];
% end 
%     save([pathstr '/' imNameC(1:end-1) '.mat']);
    load([pathstr '/' imNameC(1:end-1) '.mat']);
        

        

K0 = [1792.4519 0 1062.0368; 0 1792.4519 746.9785; 0 0 1];

[error, IAC, K] = calib_cplsvxs(cps, lss, vxs, K0);

RR = getOrientation(polyBoundaryVecCell, K, lss, lhs); 
% disp('RR1') 
% RR1{1}
% RR1{2}
% RR2 = getOrientation2(polyBoundaryVecCell, K, lss, lhs, vxs);
% disp('RR2') 
% RR2{1}
% RR2{2}
getPrjMat(imgFilename, RR, angles, K);


