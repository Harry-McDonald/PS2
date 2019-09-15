%% Extracting edges and filtering non-card objects
close all
clear all
orig = imread('C:\Users\Harry\Documents\MATLAB\YEAR 3\SEM 2\METR4202\PS2_Images\PS2 Images\Basic\Simple13.png'); %import image
im1 = rgb2gray(orig); %make grayscale
levels = multithresh(im1,5);
thresh = double(min(levels))/256;
im1 = edge(im1,'canny',thresh(1)); %find edges based on otsu threshold algorigthm
se = strel('disk',5,0);
im2 = imdilate(im1,se); %Dilate image to fill gaps based on structural element se
im3 = imfill(im2,'holes'); % fill in card face
im3 = imerode(im3,se); % erode back to normal size
im4 = edge(im3,'canny',thresh(1));
figure; imshow(im4) % This is what edge detection returns
%We can also use bwboundaries to find the coordinates of the boundaries and
%plot them over original image. This gives nicer filled in boundaries.
figure; imshow(orig)
hold on;
[B,L,n,A] = bwboundaries(im3);
%NOTE: n = number of objects
% Plotting each boundary
for ii = 1:n
    data = B{ii};
    plot(data(:,2),data(:,1),'green','LineWidth',2)
end
%% Task 2
%Next step is to filter out the obect edges
%Firstly I found the approximate area of each card using regionpropfilt and set an area range
% to be fed into bwareafilt. This worked but if the zoom was changed the
% area range would be invalid. Code commented out below.
% area_range = [62000,75000];% From the angle and distance of the photos provided this was the min and max range of the cards
% im3 = bwareafilt(im3,area_range); % Filter out all objects with area < or > area_range
%Instead I found the ratio of the minorAxisLength and majorAxisLength for
%each image found using imregionprops. If this matched the approximate
%aspect ratio of a playing card (87/56) I grabbed its area filtered out all 
%other images that werent within 5% of this area.
[im_out,properties,goodIm_index] = filterRegions(im3); %custom function made
for k = 1:length(goodIm_index)
    area = properties(goodIm_index(k)).Area;
    area_up = area+area*0.05;
    area_low = area-area*0.05;
    filt_im = bwareafilt(im3,[area_low,area_up]);
end
%Plot boundary pixels over original image. 
figure; imshow(orig)
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
%Plot centroids
for k = 1:n
    pos = centroids(k).Centroid; %2x1 coordinate vector
    plot(pos(1,1),pos(1,2),'r+')
end
for ii = 1:n
    x = centroids(ii).Centroid(1,1);
    y = centroids(ii).Centroid(1,2);
    theta = orients(ii).Orientation;
    num = num2str(ii);
    t = matlab.lang.makeValidName(num,'Prefix','t_');
    poses.(t) = [x;y;0];
    orient = matlab.lang.makeValidName(num,'Prefix','theta_');
    poses.(orient) = theta;
    N = matlab.lang.makeValidName(num,'Prefix','T_');
    poses.(N) = [cos(theta) sin(theta) 0 x;
        -sin(theta) cos(theta) 0 y;
        0 0 1 0;
        0 0 0 1];
end
%Build 2D transformation matrix based on centroid position and orientation
%The origin is from the bottom left most pixel


% for ii = 1:length(centroids)
%     pos = filt_props(ii).Centroid; %2x1 coordinate vector
%     plot(pos(1,1),pos(1,2),'r+')
% end
%% Homography
% To find the smallest distance between cards we use the boundary
% coordinates of each image (this is B from bwboundaries) and find the distance between
% each each using sqrt(x^2+y^2) we then return this distance.

% Build a for loop that iterates through all cominations of boundaries other than b1=b1. I.e
% for 3 images with boundaries it will look at (b1,b2) (b1,b3) (b2,b1)
% (b2,b3)..
num_boundaries = size(B, 1);
for b1 = 1:num_boundaries
    for b2 = 1:num_boundaries
        if b1==b2
            continue
        end
        boundary1 = B{b1};
        boundary2 = B{b2};
        boundary1x = boundary1(:,2);
        boundary1y = boundary1(:,1);
        boundary2x = boundary2(:,2);
        boundary2y = boundary2(:,1);
        for j = 1:length(boundary2y)
            %Here we are going to find the smallest x and y distance between 2 boundaries 
            bound2x = boundary2(j,2);
            bound2y = boundary2(j,1);
            all_xdist = boundary1x - bound2x;
            all_ydist = boundary1y-bound2y;
            distances = sqrt(all_xdist.^2 + all_ydist.^2);
            [min_dist, index] = min(distances);
            x1 = boundary1x(index);
            y1 = boundary1y(index);
            x2 = boundary2x(index);
            y2 = boundary2y(index);
            images = append(num2str(b1),'_to_',num2str(b2));
            minDistance = matlab.lang.makeValidName(images,'Prefix','min_dist_');
            minDistances.(minDistance) = min_dist;
        end
    end
end
% message = sprintf('Found %d boundaries', numberOfBoundaries);
% uiwait(helpdlg(message));
% Find minimum distance between each pair of boundaries
% for b1 = 1 : num_boundaries
% 	for b2 = 1 : num_boundaries
% 		if b1 == b2
% 			% Can't find distance between the region and itself
% 			continue;
% 		end
% 		bound1 = B{b1};
% 		bound2 = B{b2};
% 		bound1x = bound1(:, 2);
% 		bound1y = bound1(:, 1);
% 		x1=1;
% 		y1=1;
% 		x2=1;
% 		y2=1;
% 		overallMinDistance = inf; % Initialize.
% 		% For every point in boundary 2, find the distance to every point in boundary 1.
% 		for k = 1 : size(bound2, 1)
% 			% Pick the next point on boundary 2.
% 			bound2x = bound2(k, 2);
% 			bound2y = bound2(k, 1);
% 			% Find the closest disntance between boundary lines
% 			allDistances = sqrt((bound1x - bound2x).^2 + (bound1y - bound2y).^2);
% 			% Find closest point, min distance.
% 			[minDistance, indexOfMin] = min(allDistances);
% 			if minDistance < overallMinDistance
% 				x1 = bound1x(indexOfMin);
% 				y1 = bound1y(indexOfMin);
% 				x2 = bound2x;
% 				y2 = bound2y;
% 				overallMinDistance = minDistance(k);
% 			end
% 		end
% 		% Find the overall min distance
% 		minDistance = min(minDistance);
% 		% Report to command window.
% 		fprintf('The minimum distance from region %d to region %d is %.3f pixels\n', b1, b2, minDistance);
% 
% 		% Draw a line between point 1 and 2
% 		line([x1, x2], [y1, y2], 'Color', 'y', 'LineWidth', 3);
% 	end
% end


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
