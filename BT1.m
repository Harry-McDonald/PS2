%% Extracting edges
clear all
im1 = imread('C:\Users\Harry\Documents\MATLAB\YEAR 3\SEM 2\METR4202\PS2_Images\PS2 Images\Basic\Simple13.png'); %import image
imshow(im1)
im1 = rgb2gray(im1); %make grayscale
figure(1)
imhist(im1)
% [counts, x] = imhist(im1);
levels = multithresh(im1,5);
T = double(levels)/256;
im1 = edge(im1,'canny',T(1)); %find edges based on otsu threshold algorigthm
figure(2)
imshow(im1)
se = strel('disk',5,0);
im2 = imdilate(im1,se); %Dilate image to fill gaps based on structural element se
figure(3)
imshow(im2)
im3 = imfill(im2,'holes'); % fill in card face
im3 = imerode(im3,se); % erode back to normal size
figure(4)
imshow(im3)

B = bwtraceboundary(im3,)
