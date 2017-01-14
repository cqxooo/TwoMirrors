function img = autoCrop(imgFilename)

% Objective: automatically crop the object and its reflections in front of
%            a mirror used for camera calibration
% Author: lampson
% Date : 2015-9-24

black = [0 0 0];
red = [255 0 0];
green = [0 255 0];


I = imread(imgFilename);
[w,h,dim] = size(I);

figure(1);
imshow(I);


result = zeros(5,1);
% get 5 object color and bg color && object position

target_color = zeros(5,3);
target_pos = zeros(5,2);

for i = 1:5
    h1 = imrect;
    region = getPosition(h1);
    
    c_img = imcrop(I,region);
    [cw,ch,dim] = size(c_img);
    
    regionImg = reshape(c_img,cw*ch,dim);
    target_color(i,:) = mean(regionImg);
    
    target_pos(i,:) = [uint16(region(1) + region(3)/2) uint16(region(2))];
    
end
%object and bg color
object_color = mean(target_color(1:5,:));
%bg_color = target_color(6,:);

close(figure(1));


% reshape img as a [w*h,3] matrix
rgbM = reshape(I,w*h,dim);

% label the image into regions
gray = rgb2gray(I);
thresh = graythresh(I);
I2 = im2bw(I,thresh);

[L,regionNum] = bwlabel(I2,8);
stats = regionprops(L,'Area');
numPixels = cat(1,stats.Area);

[regionNumPixels,regionID] = sort(numPixels,'descend');

% select the biggest 20 regions
range = 20;
ave_color = zeros(range,3);
region_pos = zeros(range,2);
region_area = zeros(range,1);


% get the region belong to object
for j = 1:5
        pos = target_pos(j,:);
        target_index = pos(1)*w + pos(2);
        
        %find the region
        for order = 1:range
            currentLabel = regionID(order);

            index = find( L == currentLabel );

            if any(index == target_index)
                result(j) = currentLabel;
                break;
            end 
        end
end



img = zeros(w*h,dim);
for i = 0:regionNum
    if i == 0
%         index = find( L == 0 );
%         [len,~] = size(index);
%         img(index,:) = rgbM(index,:);
    elseif any(result == i)
        index = find( L == i );
        [len,~] = size(index);
        
        img(index,:) = rgbM(index,:);
        
    end
    

end

img = reshape(img,w,h,dim);
img = uint8(img);
gimg = double(img(:,:,1));
E = canny(gimg);
imshow(E)



