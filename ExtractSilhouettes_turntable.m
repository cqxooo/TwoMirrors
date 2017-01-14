clear all
imgFilename = {};
for i=1:34
    eval(['imgFilename{end+1} =''man/man' int2str(i)  '.jpg'';']);
end

for i=1:length(imgFilename)
    [pathstr,imNameC,imext] = fileparts(imgFilename{i});
    filename = fullfile(pathstr,[imNameC '.png']);
    im = imread(filename);
    im = double(im(:,:,1))/255;  
    cc = SelectObjectMex( im, 1 );
    polyBoundaryCell = GetBoundaryMex(cc)';  
%     ShowPoly( polyBoundaryCell, 'FaceColor', 'none', 'EdgeColor', MyPalette(i), 'LineWidth', 2);
%     axis ij,    axis equal
    filename = fullfile(pathstr,[imNameC '.sa']);
    Cp = fitBSpline(polyBoundaryCell', 10);
    writeBSplines(filename,Cp,1);
end


