function T = normalize(polyBoundaryVecCell)

%    T = normalize(Bspline)
%    T = transformation induced by data normalization
%    Bspline = input bspline points
%
%
% imgFilename=strrep(imgFilename,'png','jpg');
fig=figure(2);
x=[];
for viewLoop = 1:5,  
    hold on
    ShowPoly( polyBoundaryVecCell{viewLoop},...
        'FaceColor', MyPalette(viewLoop),  'EdgeColor',  MyPalette(viewLoop),      'FaceAlpha', 0.15);  
    x=[x;polyBoundaryVecCell{viewLoop}'];
    axis ij,    axis equal 
end
% saveas(fig,imgFilename,'jpg');
close(fig);
len = size(x,1);
mean_p = mean(x);

y(:,1) = x(:,1)-mean_p(1);
y(:,2) = x(:,2)-mean_p(2);

std_p = sqrt(sum((y(:,1)).^2+(y(:,2)).^2)/(len-1));

T = [1/std_p     0   -mean_p(1)/std_p
         0   1/std_p -mean_p(2)/std_p
         0       0           1       ];


