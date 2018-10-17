function roi = vector_roi(d, xmin, xmax, ymin, ymax)

%delete a particle if it is outside the roi at any of the time points, and
%shift the field to put the minimum at zero


time = length(d);
num_beads=length(d(1).r);

to_delete = [];

for i = 1:time %loop through times
    for k = 1:num_beads %loop through particles
        %determine if it is outside of the field, and tag for deletion
        if d(i).r(k,1)> xmax ||...
                d(i).r(k,1)< xmin ||...
                d(i).r(k,2)> ymax ||...
                d(i).r(k,2)< ymin
            to_delete = [to_delete;k];
        end
    end
end

%remove duplicates
to_delete = unique(to_delete);

%delete particles outside of the roi
for i = length(to_delete):-1:1
    for j =1:time
        d(j).r(to_delete(i),:) = [];
        d(j).dr(to_delete(i),:) = [];
    end
end

%shift the roi to put the minimum at zero

for i = 1:time
    for j = 1:size(d(1).r,1)
        d(i).r(j,1) = d(i).r(j,1) - xmin;
        d(i).r(j,2) = d(i).r(j,2) - ymin;
    end
end

roi = d;
