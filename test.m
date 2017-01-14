clc; clear all; close all;
% rgb = imread('model/baymax/baymax1.jpg');
rgb = imread('model/elephant/elephant1.jpg');
% rgb = imread('model/vase/vase1.jpg');
if ndims(rgb) == 3
    I = rgb2gray(rgb);
else
    I = rgb;
end
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

se = strel('disk', 20);
Io = imopen(I, se);

Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);

Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
fgm = imregionalmax(Iobrcbr);

se2 = strel(ones(5,5));
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);
fgm4 = bwareaopen(fgm3, 20);
It1 = rgb(:, :, 1);
It2 = rgb(:, :, 2);
It3 = rgb(:, :, 3);
It1(fgm4) = 255; It2(fgm4) = 0; It3(fgm4) = 0;
I3 = cat(3, It1, It2, It3);
bw = im2bw(Iobrcbr, graythresh(Iobrcbr));

D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;
gradmag2 = imimposemin(gradmag, bgm | fgm4);
figure('units', 'normalized', 'position', [0 0 1 1]);
subplot(2, 2, 1); imshow(bgm, []); title('分水岭变换脊线图');
subplot(2, 2, 2); imshow(fgm4, []); title('前景标记');
subplot(2, 2, 3); imshow(gradmag, []); title('梯度幅值图像');
subplot(2, 2, 4); imshow(gradmag2, []); title('修改梯度幅值图像');

L = watershed(gradmag2);
It1 = rgb(:, :, 1);
It2 = rgb(:, :, 2);
It3 = rgb(:, :, 3);
fgm5 = imdilate(L == 0, ones(3, 3)) | bgm | fgm4;
It1(fgm5) = 255; It2(fgm5) = 0; It3(fgm5) = 0;
I4 = cat(3, It1, It2, It3);
figure('units', 'normalized', 'position', [0 0 1 1]);
subplot(1, 2, 1); imshow(rgb, []); title('原图像');
subplot(1, 2, 2); imshow(I4, []); title('标记和对象边缘叠加到原图像');

Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure('units', 'normalized', 'position', [0 0 1 1]);
subplot(1, 2, 1); imshow(rgb, []); title('原图像');
subplot(1, 2, 2); imshow(Lrgb); title('彩色分水岭标记矩阵');