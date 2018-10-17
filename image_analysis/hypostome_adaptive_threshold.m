% Script to analyze opening and closing of hydra hypostome openings
% using adaptive thresholding. Input is an image sequencing with at least one opening.
% Gets user input to find the initial position of the hypostome opening.
% Created by Callen Hyland, January 2015

clear

% Enter path to folder with images
folder = './Draq5-opening-600ms/';
file_list = dir([folder '*.tif']);

% pixel-micron conversion factor
binning = 2; %binning
pix = 100/(154/binning); %microns/pixel at 10x on confocal
interval = .45; %seconds between frames

% Gaussian kernel
filt = fspecial('gaussian',[3 3], 0.5);

% size theshold for registering opening of the hypostome
cell_diam = 15; %microns
cell_area = pi*((cell_diam/2)^2); %area in microns
sz_th = cell_area/4; %minimum area for the hypostome opening to be counted
th_rng = 0.01;

% small structuring element
strel = ones(3,3);
% cell-size structuring element
cell_pix = round(cell_diam/pix);
cell_strel = strel(cell_pix, cell_pix);

%% intensity threshold for initial opening
% outline the initial opening with polygon tool

num = 100; % first frame that has mouth opening
im = im2double(imread(strcat(folder,'\',file_list(num).name)));

% select the center of hypostome opening
im_filt = imfilter(im, filt, 'replicate');
init_bw = roipoly(im_filt);
in_mask = imerode(init_bw, strel);
out_mask = imdilate(init_bw, strel);

av_out = sum(sum(((out_mask-init_bw).*im_filt)))/sum(sum(out_mask-init_bw));
av_in = sum(sum((init_bw-in_mask).*im_filt))/sum(sum(init_bw-in_mask));

init_int_th = av_in;

init_stats = regionprops(init_bw, 'Centroid', 'Area');
hyp_x = init_stats(1).Centroid(1);
hyp_y = init_stats(1).Centroid(2);

%%

hypo_area = zeros(length(file_list),1);
init_frame = 1;

% Loop through images in sequence
for i = init_frame:length(file_list)
    
    im = im2double(imread(strcat(folder,'\',file_list(i).name)));
    
    % smooth the image
    im_filt = imfilter(im, filt, 'replicate');
    
    % choose threshold for the region, search on either side of initial threshold
    th_srch = init_int_th-th_rng : 0.001 : init_int_th+th_rng;
    if i >init_frame
        diffs = zeros(size(th_srch));
        for k = 1:length(th_srch)
            test_bw = 1.-im2bw(im_filt,th_srch(k));
            test_stats = regionprops(test_bw,'Centroid');
            test_cent = [test_stats.Centroid];
            if length(test_stats)==0
                diffs(k) = 0;
            else
                dist_test = [];
                for m = 1:2:length(test_stats)*2
                    %
                    dist_test = [dist_test,sqrt((test_cent(m)-hyp_x)^2 + (test_cent(m+1)-hyp_y)^2)];
                end
                [min_test, test_ind] = min(dist_test);
                cent_min_dist = test_stats(test_ind).Centroid;
                mask = bwselect(test_bw,cent_min_dist(1),cent_min_dist(2));
                in_mask = imerode(init_bw, strel);
                out_mask = imdilate(init_bw, strel);
                av_out = sum(sum(((out_mask-init_bw).*im_filt)))/sum(sum(out_mask-init_bw));
                av_in = sum(sum((init_bw-in_mask).*im_filt))/sum(sum(init_bw-in_mask));
                diffs(k) = av_out-av_in;
            end
            
            [mx_diff, diff_ind] = max(diffs);
            
            if mx_diff > 0
                int_th = th_srch(diff_ind);
            end
            test_stats = []; test_cent = [];
        end
    else
        int_th = init_int_th;
    end
    
    %threshold and smooth edges of region
    b = im2bw(im_filt,int_th);
    b_inv = 1.-b;
    b_fill = imfill(b_inv,'holes');
    bc = imopen(b_fill,strel);
    bc_dil = imdilate(bc,strel);
    %ar = sum(sum(bc_fill))*(pix^2);
    
    cc = bwconncomp(bc_dil);
    stats = regionprops(cc, 'Area','Centroid');
    
    %find the region with the center is closest to the selected center of
    %the hypostome
    
   all_cent = [];
    if length(stats)>0
        all_cent = [stats.Centroid];
        
        dist_hyp = [];
        for j = 1:2:length(stats)*2
            dist_hyp = [dist_hyp,sqrt((all_cent(j)-hyp_x)^2 + (all_cent(j+1)-hyp_y)^2)];
        end
        [mn, mn_ind] = min(dist_hyp);
    end
    
%     all_area = [stats.Area];
%     [mx,mx_ind] = max(all_area);
    
    if length(stats)==0
        ar = 0;
    else
        %cent_max = stats(mx_ind).Centroid;
        cent_min_dist = stats(mn_ind).Centroid;
        bc_dil = bwselect(bc_dil,cent_min_dist(1),cent_min_dist(2));
        
        %remove "pockets" (uncomment this if little parts of cells at the
        %edgo of the opening are being picked up.)
        if sum(sum(bc_dil)*(pix^2)) > (6*cell_area)
            bc_dil = imerode(bc_dil, cell_strel);
            bc_dil = imdilate(bc_dil, cell_strel);
        end
        
        hypo_stats(i) = regionprops(bc_dil,'Area','Perimeter','Image');
        ar = sum(sum(bc_dil))*(pix^2);
    end
    
    border = zeros(size(im));
    if ar > sz_th
        hypo_area(i) = ar;
        bc_erd = imerode(bc_dil,strel);
        border = bc_dil-bc_erd;
        hyp_x = cent_min_dist(1);
        hyp_y = cent_min_dist(2);
    else
        hypo_area(i) = NaN;
    end
    
    imshow(imadjust(im)+border)
    title(num2str(i));
    pause(0.1)
		
		%% Save image sequence with outline of hypostome opening
		%new_dir = strcat(cd,'\','hypo_trace_movie');
		%imwrite(imadjust(im)+border,strcat(new_dir,'\',file_list(i).name),'tif');

end


%% save workspace
%save 6a-hypo-450ms
