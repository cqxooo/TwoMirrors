function polyBoundaryVecCell=loadFromTxt(imfilename, num)
[pathstr,imNameC,imext] = fileparts(imfilename);
polyBoundaryVecCell = {};
for i=1:num
    filename = fullfile(pathstr,[imNameC int2str(i) '.sa']);
    sp = loadBSplines(filename);
    fig = figure(1);
    [H, X] = drawSplines(sp(2:end,[2 1]), 10, 'g-', 1);
    polyBoundaryVecCell{i} = X(:,1:2)';
end
close(fig)