function [lineUp,lineBot,pointlUp,pointrUp,pointlBot,pointrBot] = getBitangents(O1, O2)
a = [O1'; O2'];
N = size(a,1);
length = size(O1,2);
idx = convhull(a(:,1),a(:,2));
% plot(x(idx),y(idx),'r-',x,y,'b+')
num= size(idx,1);
for i = 1:num
    if idx(i) <= length && idx(cycle(i+1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i+1,num)),:))';
        pointlUp = a(idx(i),:)';
        pointlUp = pointlUp/pointlUp(3);
        pointrUp = a(idx(cycle(i+1,num)),:)';
        pointrUp = pointrUp/pointrUp(3);
        lineUp = -sign(l(3))*l/norm(l(1:2));        
    end
    if idx(i) <= length && idx(cycle(i-1,num)) > length
        l = cross(a(idx(i),:), a(idx(cycle(i-1,num)),:))';
        pointlBot = a(idx(i),:)';
        pointlBot = pointlBot/pointlBot(3);
        pointrBot = a(idx(cycle(i-1,num)),:)';
        pointrBot = pointrBot/pointrBot(3);
        lineBot = -sign(l(3))*l/norm(l(1:2)); 
    end
end