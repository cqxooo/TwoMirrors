function  SaveAsBspline(polyBoundaryVecCell, imfilename)
[pathstr,imNameC,imext] = fileparts(imfilename); 
for i=1:length(polyBoundaryVecCell)
    filename = fullfile(pathstr,[imNameC int2str(i) '.sa']);
    Cp = fitBSpline(polyBoundaryVecCell{i}', 10);
    writeBSplines(filename,Cp,1);
end
