function [error, IAC, K] = calib_cplsvxs(cps, lss, vxs, K0)
% Calibrate with the constraints given by circular points and the rotaion
% axis, the vanishing point along the x-axis.


% w = [w1 0 w3;
%      0 w2 w4;
%      w3 w4 w5];
coef = [];
for i = 1:size(cps, 2)
    cp = cps(:,i);
    coef = [coef;     cp(1)*cp(1)+cp(2)*cp(2)   2*cp(1)       2*cp(2)                    1];     % iT w i = 0
end

for i = 1:size(lss,2)
    ls = lss(:,i);
    vx = vxs(:,i);
    coef = [ coef;
        0-ls(3)*vx(2)          ls(2)*vx(1)                ls(2)*vx(2)-ls(3)*vx(3)           ls(2)*vx(3); % ls*(w vx) = 0
        ls(3)*vx(1)            ls(3)*vx(3)-ls(1)*vx(1)          -ls(1)*vx(2)                -ls(1)*vx(3);
        -ls(2)*vx(1)+ls(1)*vx(2)    -ls(2)*vx(3)                     ls(1)*vx(3)                       0
];
end
% normalize each row
[row, column] = size(coef);

for i = 1:row
    coef(i,:)=coef(i,:)./sqrt(sum(coef(i,:).*coef(i,:)));
end

[u, d, v0] = svd(coef);

% u0=-V(2,4)/V(1,4);
% v0=-V(3,4)/V(1,4);
% s=V(4,4)+u0*V(2,4)+v0*V(3,4);
% f=sqrt(s/V(1,4));
% K=[f   0  u0
%     0  f  v0
%     0  0   1];

% v=v0;
v = sign(real(v0)).*abs(v0);
C = [v(1,end)  0   v(2,end)
    0      v(1,end)   v(3,end)
    v(2,end) v(3,end)   v(4,end)];
C = C./ abs(C(3,3));
invC = inv(C);
[tempK, p] = chol(C);
if p ~= 0
    C = -1.*C;
    [tempK, p] = chol(C);
end

K = inv(tempK);
K = K ./K(3,3)
% C is the imaged absolute conic
 IAC = C;

pp = [K0(1,3) K0(2,3)]';
f = K0(2,2);
alpha = K0(1,1)/K0(2,2);

ppN=K(1:2,3);
fN=K(2,2);
alphaN=K(1,1)/K(2,2);
epp=distPnts(pp,ppN)/sqrt(pp(1)*pp(1)+pp(2)*pp(2));
ef=abs(fN-f)/f;
eu = abs(K(1,3)-K0(1,3))/K0(1,3);
ev = abs(K(2,3)-K0(2,3))/K0(2,3);

ea=abs(alphaN-alpha)/alpha;

error= 100*[ef eu ev]


return;