clc
clear all
for i=1:21
    load(['syn' int2str(i) '.mat']);
end
K = 0; theta = 0;
K0 = [1500 0 320; 0 1500 240; 0 0 1];
theta0 = 72*pi/180;
for sigma = 0:10/20:10
    sigma
    i = 2*sigma+1;
    eval(['K = K' int2str(i) ';']);
    eval(['theta = thetas' int2str(i) ';']);
    ef = abs(K(1,1)-K0(1,1))/K0(1,1)*100
    eu = abs(K(1,3)-K0(1,3))/K0(1,3)*100
    ev = abs(K(2,3)-K0(2,3))/K0(2,3)*100
    ea = abs(theta-theta0)
end
