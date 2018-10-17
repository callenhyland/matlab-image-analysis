%Script to track nuclei in opening and closing hypostomes using manual
%input from the first image
%Last Modified: Callen Hyland, Nov. 17, 2014

clear

%folder = './Draq5-opening-600ms/';
folder = './20141112 5a-450ms/';
file_list = dir([folder '*.tif']);

new_dir = strcat(cd,'\','nuclei_track_movie');

pix = 100/77; %microns/pixel

bpass_min = 1;
bpass_max = 7;

%assumptions about the size and shape of nuclei
max_ecc = 0.95; %maximum eccentricity
ar_min = 30;    %min and max area of nuclei
ar_max = 50;
sld = 0.95;  %solidity

% Parameters for tracking
param.mem = 0; % memory; how many frames between disappearing and reappearing
param.dim = 2; % dimensions, don't change
param.good = 2; % in how many frames a particle must appear to be counted
param.quiet = 0; % suppress messages from track.m, don't change
maxdisp = 5;   % maximum displacement between frames

reg_sz = 15; %region size for looking for nuclei
SE = strel('disk', 1);
filt = fspecial('gaussian',[3 3], 0.5);

%% Show the first frame and click on all the nuclei that you want to track
% Gather an unlimited number of points until you press the return key

if ~exist('hyp_x','var')
    im = im2double(imread(strcat(folder,'\',file_list(1).name)));
    imshow(imadjust(im))
    title('select centers of nuclei. press enter when done')
    [hyp_x,hyp_y] = ginput;
else
    imshow(imadjust(im))
    hold on
    plot(hyp_x, hyp_y, 'r*')
    hold off
    title('select centers of nuclei. press enter when done')
    x = []; y = [];
    [x,y] = ginput;
    hyp_x = [hyp_x;x];
    hyp_y = [hyp_y;y];
end

% Show selected nuclei
imshow(imadjust(im))
hold on
plot(hyp_x, hyp_y, 'r*')
hold off

for i = 1:length(file_list)
    
    im = im2double(imread(strcat(folder,'\',file_list(i).name)));
    im_filt = imfilter(im,filt);
    all_regions = zeros(size(im));
    cnt = []; os = [];
    
    for j = 1:length(hyp_x)
        
        % Create subregion
        min_x = round(hyp_x(j)-reg_sz);
        max_x = round(hyp_x(j)+reg_sz);
        min_y = round(hyp_y(j)-reg_sz);
        max_y = round(hyp_y(j)+reg_sz);
        sub_im = im_filt(min_y:max_y,min_x:max_x);
        
        % Find threhsold
        int_mx = max(max(sub_im));
        auto_th = graythresh(sub_im);
        th_step = linspace(int_mx,auto_th,10);
        for k = 1:10
            im_test = im2bw(sub_im,th_step(k));
            im_test = bwselect(im_test,reg_sz,reg_sz);
            if sum(sum(im_test)) > ar_min
                level = th_step(k);
                break
            end
            level = auto_th;
        end
        
        bw = im2bw(sub_im, level);
        bw2 = im2bw(im_filt, level);
        mask = zeros(size(im));
        mask(min_y:max_y,min_x:max_x) = 1;
        mask = mask.*bw2;
%         mask = imerode(mask,SE);
%         mask = imdilate(mask,SE);
        % Select the region that is closest to the selected point
        mask = bwselect(mask,hyp_x(j),hyp_y(j));
        
        %find subpixel centroid and add to list
        os = orientation(im_filt,mask);
        tracks(j).x(i) = os.x;
        tracks(j).y(i) = os.y;
        cnt(j,1) = os.x; cnt(j,2) = os.y;
        
        all_regions = all_regions + mask;
    end
    
    % Uuse current centroids as seed points for next step
    hyp_x = cnt(:,1);
    hyp_y = cnt(:,2);
    
    imshow(all_regions), title(num2str(i))
    pause(0.1)
    %imwrite(bwperim(all_regions).*255,strcat(new_dir,'\',file_list(i).name),'tif');
    %imwrite(all_regions.*255,strcat(new_dir,'\',file_list(i).name),'tif');
    
end

%% plot tracks

for i = 1:length(tracks)
    scatter((tracks(i).x).*pix,0-(tracks(i).y).*pix,'.'), hold on
end
xlim([0 size(im,2).*pix]);
ylim([-size(im,1).*pix 0]);
xlabel('distance (microns)')
ylabel('distance (microns)')
axis equal


%% Distance traveled as function of distance from the center of hypostome
im = im2double(imread(strcat(folder,'\',file_list(1).name)));
imshow(imadjust(im))
[cen_x,cen_y] = ginput(1);
