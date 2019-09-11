%% Extracting edges
clear all
close all
orig = imread('C:\Users\Harry\Documents\MATLAB\YEAR 3\SEM 2\METR4202\PS2_Images\PS2 Images\Basic\Simple19.png'); %import image
im1 = rgb2gray(orig); %make grayscale
figure(1)
levels = multithresh(im1,5);
T = double(min(levels))/256;
im1 = edge(im1,'canny',T(1)); %find edges based on otsu threshold algorigthm
se = strel('disk',5,0);
im2 = imdilate(im1,se); %Dilate image to fill gaps based on structural element se
im3 = imfill(im2,'holes'); % fill in card face
im3 = imerode(im3,se); % erode back to normal size
figure; imshow(im3)
imageRegionAnalyzer(im3)
area_range = [62000,75000];% From the angle and distance of the photos provided this was the min and max range of the cards
im3 = bwareafilt(im3,area_range); % Filter out all objects with area < or > area_range
figure(4)
imshow(orig);
hold on;
[B,L,n,A] = bwboundaries(im3);
for ii = 1:length(B)
    data = B{ii};
    plot(data(:,2),data(:,1),'green','LineWidth',2)
end


% final_im = imoverlay(orig,card_edges,'green');
% figure(4)
% imshow(final_im)

%[im3, properties] = filterRegions(im3); %Get image properties i.e. Area of objects,axis lengths etc.


%%
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
