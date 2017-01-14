function p = intersection(L)
p = [];
for i = 1:size(L,2)-1
    for j = i+1:size(L,2)
        p = [p cross(L(:,i),L(:,j))];
    end
end
p = [p(1,:)./p(3,:);p(2,:)./p(3,:);p(3,:)./p(3,:)];