%Sholl analysis script for binary neuron image
%Callen Hyland

clear

%conversion factor
pix = 1;

%Sholl parameters
first_radius = 30;
inc = 10;
maxR = 100*inc;

path0 = './pll/thresh/';
filesi = dir([path0 '*.tif']);

%invert image
%im = 1.-(im2double(im));

%get origin for center of concentric circles
% imshow(im)
% datacursormode on
%%

for i = 1:18  %size(filesi,1)

    %stats.filename{i} = filesi(i).name;
    im = imread([path0 filesi(i).name]);
    im = 1.-(im2double(im));
    imshow(im)
    [BW, x, y] = roipoly;
    cx = round(x(1)); cy = round(y(1));
    
    im = padarray(im, [maxR maxR],1);
    im = 1.-((im));

    %column vector to hold intersection numbers
    N = maxR/inc;
    int = zeros(N,1);
    
    for j = 1:N
        
        %find points of circle with appropriate center and radius
        pts=circlepoints(cy+maxR, cx+maxR, ((j-1)*inc)+first_radius);
        %extract values of image at those points
        
        for k = 1:length(pts)
            line(k,1) = im(pts(k,1),pts(k,2));
        end
        d = diff(line);
        int(j) = sum(nonzeros(d<0));
        
    end
    pll_sholl(:,i) = int;
end


%%

h = gcf;
dcm_obj = datacursormode(h);
is = getCursorInfo(dcm_obj);
cx = is.Position(1);
cy = is.Position(2);

%pad image to insure that concentric circles do not go outside image
im = padarray(im, [maxR maxR]);

%column vector to hold intersection numbers
N = maxR/inc;
int = zeros(N,1);

for i = 1:N
    
    %find points of circle with appropriate center and radius
    pts=circlepoints(cy+maxR, cx+maxR, i*inc);
    %extract values of image at those points
    
    for j = 1:length(pts)
        line(j,1) = im(pts(j,1),pts(j,2));
    end
    d = diff(line);
    int(i) = sum(nonzeros(d<0));
    
end

[x,y] = size(blebb_sholl_20uM);
blebb_sholl_20uM(:,y+1) = int;

%%

for i = 1:12
plot(vehicle_thresh(:,i))
hold on
end