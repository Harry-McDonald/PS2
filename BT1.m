%% Extracting edges and filtering non-card objects
close all
clear all
orig = imread('C:\Users\Harry\Documents\MATLAB\YEAR 3\SEM 2\METR4202\PS2_Images\PS2 Images\Basic\Simple19.png'); %import image
im1 = rgb2gray(orig); %make grayscale
figure(1)
levels = multithresh(im1,5);
thresh = double(min(levels))/256;
im1 = edge(im1,'canny',thresh(1)); %find edges based on otsu threshold algorigthm
se = strel('disk',5,0);
im2 = imdilate(im1,se); %Dilate image to fill gaps based on structural element se
im3 = imfill(im2,'holes'); % fill in card face
im3 = imerode(im3,se); % erode back to normal size
figure; imshow(im3)
%imageRegionAnalyzer(im3)
[im_out,properties,goodIm_index] = filterRegions(im3); %custom function made
% area_range = [62000,75000];% From the angle and distance of the photos provided this was the min and max range of the cards
% im3 = bwareafilt(im3,area_range); % Filter out all objects with area < or > area_range
for k = 1:length(goodIm_index)
    area = properties(goodIm_index(k)).Area;
    area_up = area+area*0.05;
    area_low = area-area*0.05;
    filt_im = bwareafilt(im3,[area_low,area_up]);
end
%Plot boundary pixels over original image. 
figure(4)
imshow(orig);
hold on;
[B,L,n,A] = bwboundaries(filt_im);
%NOTE: n = number of objects
% Plotting each boundary
for ii = 1:n
    data = B{ii};
    plot(data(:,2),data(:,1),'green','LineWidth',2)
end
hold on;
%% Finding Pose
% Get centroid position from image properties
% Centroids is a nx1 struct containg the xy coordinate of the centroid of each card 
centroids = regionprops(filt_im,'Centroid');
orients = regionprops(filt_im,'Orientation');
for ii = 1:n
    x = centroids(ii).Centroid(1,1);
    y = centroids(ii).Centroid(1,2);
    theta = orients(ii).Orientation;
    num = num2str(ii);
    N = matlab.lang.makeValidName(num,'Prefix','T_');
    poses.(N) = [cos(theta) sin(theta) 0 x;
        -sin(theta) cos(theta) 0 y;
        0 0 1 0;
        0 0 0 1];
end
%Build 2D transformation matrix based on centroid position and orientation
%The origin is from the bottom left most pixel

%Plot centroids
for ii = 1:length(centroids)
    pos = centroids(ii).Centroid; %2x1 coordinate vector
    plot(pos(1,1),pos(1,2),'r+')
end
% for ii = 1:length(centroids)
%     pos = filt_props(ii).Centroid; %2x1 coordinate vector
%     plot(pos(1,1),pos(1,2),'r+')
% end
%% Homography
%% Random stuff 
% [H,T,R] = hough(BW);
% P = houghpeaks(H,4);
% lines = houghlines(BW,T,R,P);
% 
% for k = 1:length(lines)
%     xy=[lines(k).point1; lines(k).point2];
%     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% end

% im3 = bwareafilt(im3,[])
% 
% B = bwtraceboundary(im3,)
