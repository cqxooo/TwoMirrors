function getPrjMat(imgFilename, RR, angles, K)
PrjMat = {};
sigma = [-1 0 0;0 1 0;0 0 1];
% t = [0 0 1]'; 
% for imLoop=1:length(imgFilename)
%     [pathstr,imNameC,imext] = fileparts(imgFilename{imLoop});
%     R = RR{imLoop};
%     alpha = angles {imLoop}(1);
%     beta = angles {imLoop}(2);
%     PrjMat{1} = K*R*[eye(3)*sigma -t];
%     PrjMat{2} = K*R*[rot(2*alpha, 'y', 0) -t];
%     PrjMat{3} = K*R*[rot(2*pi-2*beta, 'y', 0) -t];
%     PrjMat{4} = K*R*[rot(2*alpha+2*beta, 'y', 0)*sigma -t];
%     PrjMat{5} = K*R*[rot(2*pi-2*alpha-2*beta, 'y', 0)*sigma -t];
%     for i=1:length(PrjMat)
%         filename = fullfile(pathstr,[imNameC int2str(i) '.pa']);
%         writePrjMat(filename, PrjMat{i});
%     end
% end
t = [0 0 -5]'; 
for imLoop=1:length(imgFilename)
    [pathstr,imNameC,imext] = fileparts(imgFilename{imLoop});
    R = RR{imLoop};
    alpha = angles {imLoop}(1);
    beta = angles {imLoop}(2);
    PrjMat{1} = K*R*[eye(3) -t];
    PrjMat{2} = K*R*[rot(-2*pi+2*alpha, 'y', 0)*sigma -t];
    PrjMat{3} = K*R*[rot(-2*beta, 'y', 0)*sigma -t];
    PrjMat{4} = K*R*[rot(-2*pi+2*alpha+2*beta, 'y', 0) -t];
    PrjMat{5} = K*R*[rot(-2*alpha-2*beta, 'y', 0) -t];
    for i=1:length(PrjMat)
        filename = fullfile(pathstr,[imNameC int2str(i) '.pa']);
        writePrjMat(filename, PrjMat{i});
    end
end