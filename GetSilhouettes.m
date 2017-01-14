function  GetSilhouettes(imfilename, num)
[pathstr,imNameC,imext] = fileparts(imfilename);
for i =1:num
    filename = fullfile(pathstr,[imNameC int2str(i) imext]);
    [image,b,c] = imread(filename,imext(2:end));
    grayImage = double(image(:,:,1));
    % canny edge detection
    S = 1;
    U = 30;
    L = 10;
    E = canny(grayImage, S, U, L);
    L = 3;
    T = 60;
    P = 30;
    C = chain(E, L, T, P);
    if C(1,2) > 1
        D = 5;
        P = 1;
        % [max_idx  no_chains   no_edgels      0     ] <---header of structure
        % [length  closed_flag    xmax       ymax    ] <---header of each chain
        len =  C(1,1); 
        C = [0 1 len 0;C];
        C(2,2) = 0;
    end
    X = [];
    len = C(1,1);
    X(:,1) = C(4:end,1);
    X(:,2) = C(4:end,2);
    figure(1);
    plot(X(:,1),X(:,2),'m-');
    axis ij;
    axis([0 2400 0 2400]);
    filename = fullfile(pathstr,[imNameC int2str(i) '.sa']);
    Cp = fitBSpline(X, 10);
    writeBSplines(filename,Cp,1);
end

