clear
clc
close all;
img= imread("night_view.png");
img = im2gray(img);
% 
img_f = double(img);
% img_f = 255.0./(1+exp(1).^(-2.*img_f+4));
% img_r = uint8(img_f);

BW1 = edge(img,'Canny');

r = Canny(img, 40, 70);

figure,imshow(img);
figure,imshow(BW1);
figure,imshow(r);
