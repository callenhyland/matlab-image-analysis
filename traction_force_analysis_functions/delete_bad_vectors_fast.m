
function [good_disp, varargout] = delete_bad_vectors_fast(d, win, varargin)
%[good_disp, bad_disp] = delete_bad_vectors(d, win, method, factor, interp_mode, graphics)
%DELETES ABERRANT VECTORS FROM A SAMPLED VECTOR FIELD, BY COMPARING EACH
%VECTOR TO A VECTOR INTERPOLATED FROM THE SURROUNDING VECTORS.
%
%REQUIRED INPUTS
%d: is structure of displacements that is output from track2disp.m.
%This is a Nx1 structure, where N is the number of time points. Each time
%point has two fields, the first is a two column array where the first
%colunm is x positions and the second column is y positions. The second
%field is a two column array where the first column is displacements in x
%and the second column is displacements in y. Each row corresponds to a
%different particle
%
%win: is a number specifying the size of the region of interest used to 
% interpolate the vector field about each sampled location. the region of
% interest is square and this number is the width of the square. this
%number should be large enough
%to include several other beads but not so large that it covers regions
%that vary greatly in displacement magnitude and direction. If it's not big
%enough there will not be enough points to interpolate.  use the optional
% graphical out put, described below, to help you pick a good value of win
%
%OPIONAL INPUTS
%'method' is a string that specifies the criterion to distiguish a bad
%vector.  The choices are:
%       'dif' or 'difference':  the magnitude of the vector difference
%       'mag' or 'magnitude':  difference of the magnitude
%       'angle':  the difference between the angles
% If you don't specify a method it will do both 'magnitude' and 'difference'.
%
%'factor' is a number specifying how much the vector is allowed to deviate 
%from the interpolated vector (default is 2). If you select 'angle' as the 
% method, then the factor should be expressed as a multiples of 45 degrees.
% (doesn't have to be an integer!)
%
%'interp_mode' is a string that specifies whether the method of
%interpolation is data gridding or surface fitting. surface fitting uses
%two subfunctions surfacefit and surfaceval. The options are:
%           'grid' : data gridding
%           'surf' or 'surface' : surface fitting
% If you don't specify and interpolation mode it will do data gridding
%
%'graphics' if you input anything in this location, this code will generate
%a handy figure (default is no figure). The figure shows vectors that were kept
%in black, vectors that were deleted in red, and vectors that didn't have 
%enough neighbros (3) in green.  This is a great option
% to turn on to help you determine which method, factor and win to use.
%
%OUTPUTS
%good_disp: all that vectors that were kept.  the structure has the same
%format as the output of track2disp
%bad_disp: all that vectors that were deleted.  the structure has the same
%format as the output of track2disp
%
%MODIFICATION HISTORY
% Created by Callen Hyland, Yale University, August 2010
%
%KNOWN ISSUES
% 1. Does not yet know how to handle vectors at the edge of the field of view
% 2. small good vectors near large bad vectors sometimes get deleted
% because the large bad vector dominates the interpolation.


%%

fig = 0;  %flag for displaying the output figure
method = [];
interp_mode = [];

% HANDLING OPTION INPUTS
if nargin == 6
    fig = 1;
    interp_mode = varargin{1,3};
    factor = varargin{1,2};
    method = varargin{1,1};
    
elseif nargin == 5
    interp_mode = varargin{1,3};
    factor = varargin{1,2};
    method = varargin{1,1};

elseif nargin == 4
    factor = varargin{1,2};
    method = varargin{1,1};

elseif nargin == 3
    factor = 2; %factor defaults to 2 if none is specified
    method = varargin{1,1};
else
    factor = 2;
end

%%

time = length(d);
num_beads=length(d(1).r);

to_delete = [];     %holds indices of vectors to delete
ignored = [];       %holds indices of vectors that do not have enough neighbors
win2 =(win/2);

%ignore the first time point because it is always zeros.
for i = 1:time %loop through times
    
    for j = 1:num_beads %loop through particles
        
        %construct window of specified size around the particle
        pos = d(i).r(j,:);
        x_win = [pos(1)-win2, pos(1)+win2];
        y_win = [pos(2)-win2, pos(2)+win2];
        
        %%%% Modified by YX to speed up (?)
        %find all particles within that window
%         field = [];

        idwin = find(abs(d(i).r(:,1)-pos(1))<win2 & ...
            abs(d(i).r(:,2) - pos(2))<win2 & ...
            (d(i).r(:,1)-pos(1)).^2 + (d(i).r(:,2) - pos(2)).^2 ~=0);
        field = [d(i).r(idwin,1:2),d(i).dr(idwin,1:2)];

%         for k = 1:num_beads
%            if d(i).r(k,1)> x_win(1)&&...
%                    d(i).r(k,1)< x_win(2) &&...
%                    d(i).r(k,2)> y_win(1) &&...
%                    d(i).r(k,2)< y_win(2) && k~=j
%                 field = [field;d(i).r(k,:), d(i).dr(k,:)];
%            end
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if size(field,1)>2
            
            switch lower(interp_mode)
                case{'grid'}
                    %use data gridding to interpolate within this window the vector at 'pos'
                    interp_x = griddata(field(:,1), field(:,2), field(:,3), pos(1), pos(2));
                    interp_y = griddata(field(:,1), field(:,2), field(:,4), pos(1), pos(2));
                    interp_disp = [interp_x interp_y];

                case{'surf', 'surface'}
                    %use surface fitting to find displacement at 'pos'
                    ax = surfacefit(field(:,1),field(:,2), field(:,3));
                    ay = surfacefit(field(:,1),field(:,2), field(:,4));
                    
                    interp_x = surfaceval(pos(1),pos(2),ax);
                    interp_y = surfaceval(pos(1),pos(2),ay);
                    interp_disp = [interp_x interp_y];
                    
                otherwise
                    %use data gridding if no method is specified
                    interp_x = griddata(field(:,1), field(:,2), field(:,3), pos(1), pos(2));
                    interp_y = griddata(field(:,1), field(:,2), field(:,4), pos(1), pos(2));
                    interp_disp = [interp_x interp_y];
            end
            
            %apply a criterion to determine if there is a discrepency
            real_disp = d(i).dr(j,:);
            real_angle = atan2(d(i).dr(j,1),d(i).dr(j,2))*180/pi;
            interp_angle = atan2(interp_x, interp_y)*180/pi;
            
            %compare vectors and delete some
            switch lower(method)
                case {'mag','magnitude'}
                    if norm(real_disp)> factor*norm(interp_disp)
                        to_delete = [to_delete;j];
                    end 
                
                case {'dif', 'difference'}
                    if norm(real_disp-interp_disp)> factor*norm(interp_disp)
                        to_delete = [to_delete;j];
                    end 
                    
                case {'angle'}
                    if abs(real_angle - interp_angle)>factor*45;
                        to_delete = [to_delete;j];
                    end
                otherwise %if no method is specified do both mag and diff
                    if norm(real_disp)> factor*norm(interp_disp)||...
                        norm(real_disp-interp_disp)> factor*norm(interp_disp) ||...
                        (abs(real_angle - interp_angle)>factor*20 && ...
                        norm(real_disp-interp_disp)> 1.5*norm(interp_disp));
                        to_delete = [to_delete;j];
                    end  
            end
        else
            ignored = [ignored;j];
        end
    end
    
end

deleted = size(to_delete,1);
display(strcat('deleting',{' '},num2str(deleted),'/',num2str(num_beads),{' '},'bad vectors'));

%make a structure of the vectors that did not have enough neighbors
ignored = unique(ignored);
ign = struct([]);
for i = length(ignored):-1:1
    for j = 1:time
        ign(j).r(i,:) = d(j).r(ignored(i),:);
        ign(j).dr(i,:) = d(j).dr(ignored(i),:);
    end
end

%delete the bad tracks, making a structure of deleted vectors
to_delete = unique(to_delete);
del = struct([]);
for i = length(to_delete):-1:1
    for j = 1:time
        del(j).r(i,:) = d(j).r(to_delete(i),:);
        del(j).dr(i,:) = d(j).dr(to_delete(i),:);
        
        d(j).r(to_delete(i),:) = [];
        d(j).dr(to_delete(i),:) = [];
    end
end


%print to the command line the number of vectors that have been deleted
deleted = length(to_delete);
display(strcat('deleting',{' '},num2str(deleted),'/',num2str(num_beads),{' '},'bad vectors'));

%print to the command line the number of vector that did not have enough
%neighbors
if ~isempty(ign)
    display(strcat(num2str(length(ignored)),{' '},'vectors did not have enough neighbors'));
end

%Now display the figure if the user has asked for it
if fig
    figure(100)
    for i = 1%2:time
        drf=mean(d(i).dr); scl = 10;
        quiver((d(i).r(:,1)),(d(i).r(:,2)),scl*(d(i).dr(:,1)-drf(1)),scl*(d(i).dr(:,2)-drf(2)),0,'black');
        hold on
        if ~isempty(del)
            quiver((del(i).r(:,1)),(del(i).r(:,2)),scl*(del(i).dr(:,1)-drf(1)),scl*(del(i).dr(:,2)-drf(2)),0,'red');
        end
        if ~isempty(ign)
            quiver((ign(i).r(:,1)),(ign(i).r(:,2)),scl*(ign(i).dr(:,1)-drf(1)),scl*(ign(i).dr(:,2)-drf(2)),0,'green');
        end
        hold off
    end
end

%OUTPUTS
good_disp = d;
if nargout>1
    varargout{1,1} = del;
end

function a = surfacefit(X,Y,U)

%X is a vector of x coordinates
%Y is a vector containing corresponding y coordingates
%U is a vector of values corresponding to the x and y coordinates
%all three of these vectors must be the same length

%break if there are not enough data points to solve
n = size(X,1);
if n<6
    error('not enough data points to fit surface')
end

%construct matrix P
P = zeros(n,6);
P(:,1) = 1;
P(:,2) = X;
P(:,3) = Y;
P(:,4) = X.^2;
P(:,5) = Y.^2;
P(:,6) = X.*Y;

a = pinv(P'*P)*P'*U;

function z = surfaceval(x,y,a)

%x and y are the coordinates of the point to evaluate the function at
%a is a column vector containing exactly six coefficients for the equation
%of the surface. a is the output from the 'surfacefit' function.

z = a(1) + a(2)*x + a(3)*y + a(4)*(x^2) + a(5)*(y^2) + a(6)*x*y;