%Date         : 2015-5-20
%Author       : lampson
%Description  : Mean-shift algorithm for pixel clustering in an image
%Input        : a foot image
%Output       : clustered image

function mean_shift(path)
% function mean_shift
    path = 'hand\DSC_6356.jpg';
    
    img = imread(path);
    img = imresize(img,0.5);
    
    [width height channel] = size(img);
    
    %initiate labels for groups
    
    labels = zeros(width,height) - 1;
    
    radius = 40;
    seedFlag = true;
    ite = 0;
    shift = [100 100 100];
    
    while(ismember(-1,labels))
        % the (ite) th group
        ite = ite + 1;
        img = double(img);
        
        while(norm(shift) > 20)
            
            tempShift = [0 0 0];
            count = 0;
            %iterate the whole image pixel by pixel
            for i = 1:width
                for j = 1:height
                    % check whether the current pixel is grouped
                    if labels(i,j) == -1 && seedFlag
                        seed = img(i,j,:);
                        labels(i,j) = ite;
                        seedFlag = false;
                    % shift vector
                    elseif labels(i,j) == -1
                        distance = dist(img(i,j,:),seed);
                        if distance < radius
                           labels(i,j) = ite;
                           tempShift = shif( img(i,j,:),seed ) + tempShift;
                           
                           img(i,j,:) = red;
                           
                           count = count + 1;
                        end
                    end
                end
            end
            
            shift = tempShift/count;
            seed = move(seed,shift);
            
        end
        % end of current iteration
        img = uint8(img);
        figure;
        imshow(img);
        %select new seed
        seedFlag = true;
        shift = [100 100 100];
        %break;
    end
    
end

function distance = dist(p1,p2)
    pixel_1 = [p1(1,1,1) p1(1,1,2) p1(1,1,3)];
    pixel_2 = [p2(1,1,1) p2(1,1,2) p2(1,1,3)];
    distance = norm(pixel_2-pixel_1);
end

function shift = shif(p1,p2)
    pixel_1 = [p1(1,1,1) p1(1,1,2) p1(1,1,3)];
    pixel_2 = [p2(1,1,1) p2(1,1,2) p2(1,1,3)];
    shift = pixel_1-pixel_2;
end

function seed = move(s1,shift)
    seed(1,1,1) = s1(1,1,1) + shift(1);
    seed(1,1,2) = s1(1,1,2) + shift(2);
    seed(1,1,3) = s1(1,1,3) + shift(3);
end

function red_pixel = red
    red_pixel(1,1,1) = 255;
    red_pixel(1,1,2) = 0;
    red_pixel(1,1,3) = 0;
end