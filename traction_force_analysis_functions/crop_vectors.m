function [cropped_d, varargout] = crop_vectors(d, varargin)
%function [cropped_d, ROI_X, ROI_Y] = crop_disp(d, roi_x, roi_y, method)

%REQUIRED INPUTS:
%d is structure of displacements that is output from track2disp.m.
%This is a Nx1 structure, where N is the number of time points. Each time
%point has two fields, the first is a two column array where the first
%colunm is x positions and the second column is y positions. The second
%field is a two column array where the first column is displacements in x
%and the second column is displacements in y. Each row corresponds to a
%different particle. For each of the different time points, vectors
%specifying the positions and the displacements have to be the same length!
%
%this can be used to crop vector fields from multiple time points to the
%same size and make sure there are no vectors in one timepoint that do not
%appear in all the other timepoints
%
%OPTIONAL INPUTS
%roi_x and roi_y are vectors containing the x and y coordinates of the
%corners of the roi polygon. These vectors must be the same length. If you
%do not specify the corners of the roi, the function will generate a
%figure and an interactive tool that will allow you to select a region of
%interest. coordinates of the corners of the roi can be returned as an
%output

%'method' is a string that specifies what to do with the vectors that are
%outside the roi. options are:
%    'delete' : delete all vectors outside the roi and shifts the field so
%    that the minimum x and minimum y is set to zero. This is the default.
%
%   'delete_inv' is the same as 'delete' but it deletes all of the vectors
%   INSIDE roi. There is no shift in the output for this method.
%
%    'zero' : simply sets all of the vectors outside the field of view to
%    zero.
%
%    'zero_inv' is the same as zero except all of the vectors inside the
%    roi are set to zero.

%REQUIRED OUTPUTS:
%cropped_d is a structure of the same format as the output of track2disp, either
%cropped and shifted or with all vectors outside the roi set to zero.
%
%OPTIONAL OUTPUTS:
%ROI_X and ROI_Y will return two vectors containing the x and y coordinates
%of the corners of the roi. This is useful if you have used the interactive
%tool and want to reuse the same region in a susequent analysis
%
%MODIFICATION HISTORY
%Created by Callen Hyland, October 2010
%Fixed 'zero' and 'zero_inv' methods, added autoscaling for vectors, CH 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GET INPUTS
%flag to display the interactive tool
interactive = 0;
method = [];

if nargin == 4
    if isempty(varargin{1,1}) || isempty(varargin{1,2})
        method = varargin{1,3};
        interactive = 1;
    else
        roi_x = varargin{1,1};
        roi_y = varargin{1,2};
        method = varargin{1,3};
    end
    
elseif nargin == 3
        roi_x = varargin{1,1};
        roi_y = varargin{1,2};
        method = 'delete';
else
    interactive = 1;
    method = 'delete';
end


%number of timepoints
time = length(d);

%choose time point for display: 
%1 if there is only one time point, or if first time point in non-zero
if time==1 || sum(sum(d(1).dr))~=0
    k=1;
else
    k = 2;
end

%display the interactive tool if the user has not specified a polygon
if interactive == 1
    x_pos = d(k).r(:,1);
    y_pos = d(k).r(:,2);
    x_disp = d(k).dr(:,1);
    y_disp = d(k).dr(:,2);
    x_max = max(x_pos);
    y_max = max(y_pos);
    imsize = 512;
    figure
    blank_plot = zeros(imsize,imsize);
    imshow(blank_plot); hold on
    %rescale positions to fit on blank image
    sc = imsize/x_max;
    %quiver((x_pos./x_max).*imsize, (y_pos./y_max).*imsize, x_disp, y_disp, 0, 'white')
    quiver(x_pos.*sc, y_pos.*sc, x_disp, y_disp, 1, 'white')
    hold off
    [BW, roi_x, roi_y] = roipoly;
    
    %undo scaling on polygon corners
    for i = 1:length(roi_x)
        roi_x(i) = (roi_x(i)./imsize).*x_max;
        roi_y(i) = (roi_y(i)./imsize).*y_max;
    end
end


%number of vectors at each timepoint
num_beads=length(d(1).r);
in = zeros(num_beads,time);
xmn = min(roi_x);
ymn = min(roi_y);

for i = 1:time
    in(:,i) = inpolygon(d(i).r(:,1),d(i).r(:,2), roi_x, roi_y);
end

all_times = sum(in,2);
to_delete = [];

%delete or keep depending on the user-specified method
switch lower(method)
    case {'delete'}
        %if it's not inside the polygon at all times, tag it for deletion
        for i = 1:num_beads
            if all_times(i)< time
                to_delete = [to_delete;i];
            end
        end
        %remove vectors tagged for deletion
        for i = length(to_delete):-1:1
            for j =1:time
                d(j).r(to_delete(i),:) = [];
                d(j).dr(to_delete(i),:) = [];
            end
        end
        %shift the roi to put the minimum at zero
        for i = 1:time
            for j = 1:size(d(1).r,1)
                d(i).r(j,1) = d(i).r(j,1) - xmn;
                d(i).r(j,2) = d(i).r(j,2) - ymn;
            end
        end
        
    case {'delete_inv'}
        %if it's not outside the polygon at all times tag it for deletion
        for i =1:num_beads
            if all_times(i)>0
                to_delete = [to_delete;i];
            end
        end
        %remove vectors tagged for deletion
        for i = length(to_delete):-1:1
            for j =1:time
                d(j).r(to_delete(i),:) = [];
                d(j).dr(to_delete(i),:) = [];
            end
        end

    case {'zero'}
        %if it's not inside polygon at all times set magnitude to zero
        for i = 1:num_beads
            if all_times(i)< time
                for j = 1:time
                    d(j).dr(i,1) = 0;
                    d(j).dr(i,2) = 0;
                end
            end
        end
    
    case {'zero_inv'}
        %if it's not outside polygon at all times set magnitude to zero
        for i = 1:num_beads
            if all_times(i)>0
                for j = 1:time
                    d(j).dr(i,1) = 0;
                    d(j).dr(i,2) = 0;
                end
            end
        end
        
end

%OUTPUTS
if nargout==3
    varargout{1,1} = roi_x;
    varargout{1,2} = roi_y;
end

close all

cropped_d = d;
