function [lineUp,lineBot,pointUp,pointBot,flag] = getOuterTengents(pbv,e)
% this function is try to get the two outer tangent lines. flag == 0 means 
% either there is no bottem outer tangent line or there is no top outer tangent line.
% flag == 1 means it's OK.

e=e./e(3);

x = pbv(:,1);
y = pbv(:,2);
idx = convhull(x,y);
num= size(idx,1);
up = -1;
down = -1;
numdown = -1;
numup = -1;
for i = 1:num
    p = [x(idx(i)) y(idx(i)) 1]';
    l = cross(p,e);
    r = [x(idx) y(idx) ones(num,1)]*l;
    
    if size(find(r>=0),1) > numdown
        down = idx(i);
        numdown = size(find(r>=0),1);
        pointBot = pbv(down,:)';
        lineBot = -sign(l(3))*l/norm(l(1:2));
    end
    
    if size(find(r<=0),1) > numup
        up = idx(i);
        numup = size(find(r<=0),1);
        pointUp = pbv(up,:)';
        lineUp = -sign(l(3))*l/norm(l(1:2));

    end
end
% if pointUp(2)>pointBot(2)
%     tmp = pointUp;
%     pointUp = pointBot;
%     pointBot = tmp;
%     tmp = lineUp;
%     lineUp = lineBot;
%     lineBot = tmp;
% end
flag = 1;

return;

