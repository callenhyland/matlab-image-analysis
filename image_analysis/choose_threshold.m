function thresh = choose_threshold(im, sz_range)
%function thresh = choose_threshold(im, sz_range)
%
%This function finds a threshold for an image of bright objects on a dark
%background. The threshold maximizes the number of objects within a
%user-specified size range. It is helpful when tracking particles as an
%input to bwlabel.m or pkfnd.m
%
%INPUTS
%im is the input image you would like to find a threshold for- it should be
%grayscale and can be any bit depth
%
%sz is a two element vector containing minimum and maximum values for the
%diameter of the objects you are looking for. EX: [6 13]
%
%OUTPUTS
%thresh is value of the threshold that maximizes the number of objects
%within your size range. The function will also display a binary image
%created by applying the treshold to the input image.
%
%Created by Callen Hyland, August 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%area range of objects
min_area = pi*((0.5*sz_range(1))^2);
max_area = pi*((0.5*sz_range(2))^2);

%intensity range of input image
mx = max(max(im));
mn = min(min(im));

if mx<1
   im_in = im.*1000; 
   mx = mx*1000; mn = mn*1000;
else
    im_in = im;
end

center = (mx+mn)/2;
diff = mx-center;
part = [0 0];
thresh = [0 0];

while diff>1
   thresh(1) = center + (diff*0.5);
   thresh(2) = center - (diff*0.5);
   
   for i = 1:2
       mask = im_in > thresh(i);
       r = regionprops(mask,'Area');
       n=0;
       for j = 1:length(r)
           if r(j).Area > min_area && r(j).Area < max_area
               n=n+1;
           end
       end
       part(i) = n;
   end
    
   if part(1)>part(2)
       diff = thresh(1)-center;
       center = thresh(1);
   else
       diff = center-thresh(2);
       center = thresh(2);
   end
    
end

thresh = center;

if max(max(im))<1
    thresh = thresh/1000;
end

bw = im>thresh;
%imshow(bw)

