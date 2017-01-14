clear all
imgFilename = {};
for i=1:2
    eval(['imgFilename{end+1} =''formike/baymax/baymax' int2str(i)  '.jpg'';']);
end
for i=1:2
    eval(['imgFilename{end+1} =''formike/elephant/elephant' int2str(i)  '.jpg'';']);
end
for i=1:2
    eval(['imgFilename{end+1} =''formike/girl/girl' int2str(i)  '.jpg'';']);
end

for i=1:length(imgFilename)
    [pathstr,imNameC,imext] = fileparts(imgFilename{i});
    for num = 1:5,    
        filename = fullfile(pathstr,[imNameC int2str(num) '.png']);
        im = imread(filename);
        im = double(im(:,:,1))/255;  
        cc = SelectObjectMex( im, 1 );
        polyBoundaryCell = GetBoundaryMex(cc)';  
%         ShowPoly( polyBoundaryCell, 'FaceColor', 'none', 'EdgeColor', MyPalette(num), 'LineWidth', 2);
%         axis ij,    axis equal
        filename = fullfile(pathstr,[imNameC int2str(num) '.sa']);
        Cp = fitBSpline(polyBoundaryCell', 10);
        writeBSplines(filename,Cp,1);
    end
end


