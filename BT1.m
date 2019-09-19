%% Initialise
close all
clear all
clc

path = addpath(genpath("PS2 Images"));
load cameraParamsBasic cameraParams
imCount = 20;           % The number of images you want to load
imNames{imCount} = {};  % Empty array of imCount elements
prefix = "Simple";       % The prefix to use for loading images

for n = 1:imCount
   imNames{n} = ['Simple', num2str(n), '.png'];
end
% Change image_num to the required image
response = input('Please enter the image number (1-25): ');
orig  = imread(imNames{response});
orig  = undistortImage(orig, cameraParams);

%% Task 1: Extracting edges and filtering non-card objects
fprintf('\n---------------------------------- Task 1 ----------------------------------\n')
fprintf('See figure 1 for object edge detection \nSee figure 2 for edge overlay onto original image \n')
im1 = rgb2gray(orig); %make grayscale
levels = multithresh(im1,5);
thresh = double(min(levels))/256;
im1 = edge(im1,'canny',thresh(1)); %find edges based on otsu threshold algorigthm
se = strel('disk',5,0);
im2 = imdilate(im1,se); %Dilate image to fill gaps based on structural element se
im3 = imfill(im2,'holes'); % fill in card face
im3 = imerode(im3,se); % erode back to normal size
im4 = edge(im3,'canny',thresh(1));
figure; 
% subplot(1,2,1);
imshow(im4)
title('Task 1: Edges - Including non-card objects');
hold on
% This is what edge detection returns
%We can also use bwboundaries to find the coordinates of the boundaries and
%plot them over original image. This gives nicer filled in boundaries.

figure;
imshow(orig)
hold on;
[B,L,n,A] = bwboundaries(im3);
%NOTE: n = number of objects
% Plotting each boundary
for ii = 1:n
    data = B{ii};
    plot(data(:,2),data(:,1),'green','LineWidth',2)
    title(' Task 1: Edge Overlay - Including non-card objects')
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
%aspect ratio of a playing card (87/56) I grabbed its area and filtered out all 
%other images that werent within 5% of this area.
fprintf('\n---------------------------------- Task 2 ----------------------------------\n')
fprintf('See figure 3 for card edge detection that filters out non-card objects.\n')
[im_out,properties,goodIm_index] = filterRegions(im3); %custom function made
for k = 1:length(goodIm_index)
    area = properties(goodIm_index(k)).Area;
    area_up = area+area*0.1;
    area_low = area-area*0.1;
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
title('Task 2: Edge Overlay - Non-card objects filtered out')
hold on;
%% Task 3: Finding Pose
fprintf('\n---------------------------------- Task 3 ----------------------------------\n')

[B,L,n,A] = bwboundaries(filt_im);
numberOfCards = n;
fprintf('There are %d cards in this image.\n\n',n)
% Get centroid position from image properties
% Centroids is a nx1 struct containg the xy coordinate of the centroid of each card 

Poses = regionprops('Table',filt_im,{'Centroid','Orientation'});
centroids = regionprops(filt_im,'Centroid');
orients = regionprops(filt_im,'Orientation');
%Plot centroids and origin
figure;
imshow(orig)
hold on
for k = 1:n
    x_0 = 0;
    y_0 = 0;
    pos = centroids(k).Centroid; %2x1 coordinate vector
    x = pos(1,1);
    y = pos(1,2); 
    plot(x,y,'+','MarkerSize',15, 'DisplayName', ['Card',num2str(k)])
    text(x+ 30,y+30,['(',num2str(round(x,2)),', ', num2str(round(y,2)),')'],'color','b')
end
title('Object Index')
legend;
hold off
% This for loop finds the poses of each card based on the orientation and
% centroid position
for ii = 1:n
    x = centroids(ii).Centroid(1,1);
    y = centroids(ii).Centroid(1,2);
    %Build cell of centroid positions for plotting later
    centroid_pos{ii,1} = x;
    centroid_pos{ii,2} = y;
    theta = orients(ii).Orientation;
    num = num2str(ii);
    t = matlab.lang.makeValidName(num,'Prefix','Translation_');
    translations.(t) = [x;y];
    orient = matlab.lang.makeValidName(num,'Prefix','Orientation_');
    orientations(ii).(orient) = theta;
    N = matlab.lang.makeValidName(num,'Prefix','Transform_');
    T_forms.(N) = [cos(theta) sin(theta) 0 x;
        -sin(theta) cos(theta) 0 y;
        0 0 1 0;
        0 0 0 1];
end
if n == 1
    fprintf('There is 1 card in this image.\n\n',n)
    fprintf("The translation of the card's centroid is: (in pixels)\n")
    struct2table(translations)
    fprintf('The orientation of the card, measured as the angle from the major axis to the positve x axis, is: (in degrees) \n',n)
    struct2table(orientations)
    fprintf('The transformation matrix for this card is:\n')
    struct2table(T_forms)
end
if n>1
    fprintf('There are %d cards in this image.\n\n',n)
    fprintf('The translation of the centroids for cards 1 - %d are: (in pixels)\n',n)
    struct2table(translations)
    fprintf('The orientations of cards 1 - %d are: (in degrees)\n',n)
    struct2table(orientations)
    fprintf('The transformations matrices for cards 1 - %d are:\n',n)
    struct2table(T_forms)
end
% for ii = 1:n
%     x = centroids(ii).Centroid(1,1);
%     y = centroids(ii).Centroid(1,2);
%     %Build cell of centroid positions for plotting later
%     centroid_pos{ii,1} = x;
%     centroid_pos{ii,2} = y;
%     theta = orients(ii).Orientation;
%     num = num2str(ii);
%     t = matlab.lang.makeValidName(num,'Prefix','t_');
%     poses.(t) = [x;y;0];
%     orient = matlab.lang.makeValidName(num,'Prefix','theta_');
%     poses.(orient) = theta;
%     N = matlab.lang.makeValidName(num,'Prefix','T_');
%     poses.(N) = [cos(theta) sin(theta) 0 x;
%         -sin(theta) cos(theta) 0 y;
%         0 0 1 0;
%         0 0 0 1];
% end

%% Homography
% To find the smallest distance between cards we use the boundary
% coordinates of each image (this is B from bwboundaries) and find the distance between
% each each using sqrt(x^2+y^2) we then return this distance.

% Build a for loop that iterates through all cominations of boundaries other than b1=b1. I.e
% for 3 images with boundaries it will look at (b1,b2) (b1,b3) (b2,b1)
% (b2,b3)..
figure; 
imshow(orig)
hold on
legend('Location','southeast');
% for ii = 1:n
%     data = B{ii};
%     plot(data(:,2),data(:,1),'green','LineWidth',2,'HandleVisibility','off')
% end
for j = 1:n
    plot(centroid_pos{j,1},centroid_pos{j,2},'+','DisplayName',['Card ', num2str(j)])
end
num_boundaries = size(B, 1);
combo =[];% Initialise combination list, this will be used to prevent repetitions of combinations
for b1 = 1:num_boundaries
    if num_boundaries == 1 
        minDistances = 0;
        break
    end
    for b2 = 1:num_boundaries
        cont = 0; % Initialise continue variable - this will be used trigger the next iteration of combinations if there is a repetition
        % Gather list of combinations, check their are no repetitions
        combo = [combo;[b1,b2]];
        sz_combo = size(combo);
        combo_rows = sz_combo(1);
        for c = 1:combo_rows
            if combo(c,2) == b1 && combo(c,1)==b2
                cont = 1;
                break
            end
        end
        % Check boolean of cont to see whether we need to skip this
        % iteration
        if cont == 1
            continue
        end
        
        boundary1 = B{b1};
        boundary2 = B{b2};
        boundary1x = boundary1(:,2);% All boundary 1 x values
        boundary1y = boundary1(:,1);% All boundary 1 y values
        boundary2x = boundary2(:,2);% All boundary 2 x values
        boundary2y = boundary2(:,1);% All boundary 2 y values
        overallMinDist = inf; %This is initialising the minimum distance it will be refined by the following loop
        for j = 1:length(boundary2y)
            %Here we are going to find the smallest x and y distance between 2 boundaries
            % Therefore we need to iterate through every x and y difference
            % and find which is the smallest
            bound2x = boundary2(j,2); %second boundary x value being checked
            bound2y = boundary2(j,1); % second boundary y value being checked
            all_xdist = boundary1x - bound2x; % Difference between boundary 2 x value and all boundary 1 x values.
            all_ydist = boundary1y - bound2y; % Difference between boundary 2 y value and all boundary 1 y values.
            distances = sqrt(all_xdist.^2 + all_ydist.^2); % Hypotenuse of x and y distances
            [min_dist, index] = min(distances); % Minimum distance for boundary 2 jth x and y values
            % Need to check if jth min_dist is smaller than (j-1)th min_dist 
            if min_dist < overallMinDist %If jth min_dist is < (j-1)th min_dist, enter and record as overallMinDist
                x1 = boundary1x(index);
                y1 = boundary1y(index);
                x2 = bound2x;
                y2 = bound2y;
                overallMinDist = min_dist;
            end
        end
        images = append(num2str(b1),'to',num2str(b2));
        minDistance = matlab.lang.makeValidName(images,'Prefix','minDist');
        dist_mm = pixel2mm(overallMinDist,filt_im);
        minDistances.(minDistance) = dist_mm;
        lines = append(num2str(b1),' to ',num2str(b2));
        % Draw a line between point 1 and 2
		line([x1, x2], [y1, y2],'color',rand(1,3),'LineWidth', 3,'DisplayName',[num2str(dist_mm),' mm']);
    end
end
%minDistances
title('Minimum Paths Between Cards')

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
%% FUNCTIONS
function [BW_out,properties,goodIm_index] = filterRegions(BW_in)
%filterRegions  Filter BW image using auto-generated code from imageRegionAnalyzer app.
%  [BW_OUT,PROPERTIES] = filterRegions(BW_IN) filters binary image BW_IN
%  using auto-generated code from the imageRegionAnalyzer app. BW_OUT has
%  had all of the options and filtering selections that were specified in
%  imageRegionAnalyzer applied to it. The PROPERTIES structure contains the
%  attributes of BW_out that were visible in the app.

% Auto-generated by imageRegionAnalyzer app on 11-Sep-2019
%---------------------------------------------------------

BW_out = BW_in;

% Get properties.
properties = regionprops(BW_out, {'Area', 'MajorAxisLength', 'MinorAxisLength', 'Orientation'});
aspect_list = [];
goodIm_index =[];
card_aspect = 87/56;
card_aspect_up = card_aspect + card_aspect*0.05;
card_aspect_low = card_aspect-card_aspect*0.05;
for ii = 1:length(properties)
    minor_length = properties(ii).MinorAxisLength;
    major_length = properties(ii).MajorAxisLength;
    aspect_ratio = major_length/minor_length;
%     aspect_list = [aspect_list, aspect_ratio];
    if aspect_ratio <= card_aspect_up && aspect_ratio >= card_aspect_low
        goodIm_index = [goodIm_index,ii];
    end
end
end
%Pixel to mm converter
function length_mm = pixel2mm(length_pixels,image)
    props = regionprops('table',image,'MajorAxisLength');%Card in image major axis length (pixels)
    majorCard_mm = 87; %UQ card major axis length (mm)
    avg_major_px = mean(props.MajorAxisLength);
    %majorCard_px = props.MajorAxisLength;
    mm_per_pixel = majorCard_mm/avg_major_px; %mm per pixel ratio
    length_mm = length_pixels * mm_per_pixel;
end
%Majorlength = properties(1).Area
% Uncomment the following line to return the properties in a table.
% properties = struct2table(properties);