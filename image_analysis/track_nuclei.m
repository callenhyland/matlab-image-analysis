%% Script to track nuclei in image sequences of opening and closing hypostomes
%% Created by Callen Hyland, January 2015

clear

% Path to image sequence
folder = './slide1 Draq5-DS/';
file_list = dir([folder '*.tif']);

% Assumptions about the size and shape of nuclei
max_ecc = 0.97; %maximum eccentricity
ar_min = 35;    %min and max area of nuclei
ar_max = 110;
sld = 0.9;  %solidity

% Parameters for tracking
param.mem = 0; % memory; how many frames between disappearing and reappearing
param.dim = 2; % dimensions, don't change
param.good = 2; % in how many frames a particle must appear to be counted
param.quiet = 0; % suppress messages from track.m, don't change
maxdisp = 10;   %maximum displacement between frames

% Step between thresholds
th_step = 0.0001;

% Gaussian kernel
filt = fspecial('gaussian',[3 3], 0.5);

%% Method for finding nuclei that uses a single threshold for each image

wtcentr = [];

for i = 1:length(file_list)
    
    im = im2double(imread(strcat(folder,'\',file_list(i).name)));
    im_filt = bpass(im,1,13);
    b_good = zeros(size(im_filt));
    
    th = choose_threshold(im_filt, [7 15]);    %threshold for image
    b = im2bw(im_filt,th);
    
    stats = regionprops(b, 'Eccentricity','Area','Centroid','Solidity');
    for k = 1:length(stats)
        if stats(k).Eccentricity<max_ecc...
                && ar_min<stats(k).Area && stats(k).Area<ar_max...
                && stats(k).Solidity > sld
            centr = stats(k).Centroid;
            b_good = b_good + bwselect(b,centr(1),centr(2));
        end
    end
    stats = [];
    
    im_label = bwlabel(b_good);
    
    os = [];
    for j = 1:max(max(im_label))
        mask = im_label == j;
        os = [os; orientation(im_filt,mask)];
    end
    
    c = [os.x ; os.y]';
    c(:,3) = i;
    wtcentr = [wtcentr;c];
    c = [];
    
    imshow(b_good)
    pause(0.1)
    nuclei(:,:,i) = b_good;
end

track_data = track(wtcentr, maxdisp, param); % tracks
dt=track2disp(track_data);

%% Manually select the background of the image
im = im2double(imread(strcat(im_dir,'\',file_list(2).name)));
im_filt = imfilter(im, filt, 'replicate');
imshow(imadjust(im_filt));
[mask, roi_x, roi_y] = roipoly;
bkgrd = sum(sum(mask.*im_filt))/sum(sum(mask));


%% Method for finding nuclei that uses multiple thresholds

%im_min = min(min(im_filt));
im_min = bkgrd;
im_max = max(max(im_filt));
int = im_min:th_step:im_max;

for i = 1:length(file_list)
    
    im = im2double(imread(strcat(im_dir,'\',file_list(i).name)));
    im_filt = bpass(im,1,13);
    b_good = zeros(size(im_filt));
    
    for j = 1:length(int)
        
        th = int(j);    %threshold for image
        b = im2bw(im_filt,th);
        
        stats = regionprops(b, 'Eccentricity','Area','Centroid','Solidity');
        for k = 1:length(stats)
            if stats(k).Eccentricity<max_ecc...
                    && ar_min<stats(k).Area && stats(k).Area<ar_max...
                    && stats(k).Solidity > 0.9
                centr = stats(k).Centroid;
                b_good = b_good + bwselect(b,centr(1),centr(2));
            end
        end
        stats = [];
    end
    imshow(b_good)
    nuclei(:,:,i) = b_good;
end
