function [alpha beta]=getAngle(v,lh,ls,cp, T)

v = inv(T)*v;
v = [v(1,:)./v(3,:); v(2,:)./v(3,:);v(3,:)./v(3,:)];
vz=cross(lh, ls);
vz = vz./vz(3); 
v1=v(:,1);
v2=v(:,2);
v3=v(:,3);
v4=v(:,4);

alpha = pi/2-abs(log(crossRatio(vz, v1, cp(:,1), cp(:,2)))/2*sqrt(-1));
beta = pi/2-abs(log(crossRatio(v2, vz, cp(:,1), cp(:,2)))/2*sqrt(-1));
EvalPrint('(alpha) * (180/pi)')
EvalPrint('(beta) * (180/pi)')
EvalPrint('(alpha+beta) * (180/pi)')


% alpha_ = log(crossRatio(vz, v1, cp(:,1), cp(:,2)))/2*sqrt(-1);
% alpha = pi/2-alpha_;
% beta_ = log(crossRatio(v2, vz, cp(:,1), cp(:,2)))/2*sqrt(-1);
% beta = pi/2-beta_;
% EvalPrint('(alpha_) * (180/pi)')
% EvalPrint('(beta_) * (180/pi)')
% EvalPrint('(alpha) * (180/pi)')
% EvalPrint('(beta) * (180/pi)')
% EvalPrint('(alpha+beta) * (180/pi)')


% v1=v(2,1);
% v2=v(2,2);
% v3=v(2,3);
% v4=v(2,4);
% 
% q1(2) = (v1*(v2+v3)-2*v2*v3)/(2*v1-v2-v3);
% q1(1) = -(lh(3)+lh(2)*q1(2))/lh(1);
% q1(3) = 1;
% q1 = q1';
% 
% q2(2) = (v2*(v1+v4)-2*v1*v4)/(2*v2-v1-v4);
% q2(1) = -(lh(3)+lh(2)*q2(2))/lh(1);
% q2(3) = 1;
% q2 = q2';
% 
% alpha = abs(log(crossRatio(q1, vz, cp(:,1), cp(:,2)))/2*sqrt(-1));
% beta = abs(log(crossRatio(q2, vz, cp(:,1), cp(:,2)))/2*sqrt(-1));

