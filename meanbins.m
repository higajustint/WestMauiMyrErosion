function disc_all = meanbins(binby,binof,binedges)
%MEANBINS Summary of this function goes here
%   Input
%   binby = x-axis variable to bin dependent variable
%   binof = y-axis variable to be binned based on x-axis bin edges
%   binedges = bin edges
%   Output
%   disc_all = col 1 is the midpoint of bins, col 2 is the mean y axis
%   value of bins, col 3 is the minimum y-axis value in bins, col 4 is the
%   maximum y-axis value of bins

binall = discretize(binby,binedges);
binall(:,2) = binof;

binmid = NaN(size(1,length(binedges)-1));
for i = 1:length(binedges)-1
    binmid(i) = (binedges(i) + binedges(i+1))/2;
end

disc_mean = NaN(size(1,max(binall(:,1))));
disc_min = NaN(size(1,max(binall(:,1))));
disc_max = NaN(size(1,max(binall(:,1))));
for i = 1:max(binall(:,1))
    disc_idx = find(binall(:,1) == i);
    disc_mean(i) = mean(binall(disc_idx,2));
    disc_min(i) = min(binall(disc_idx,2));
    disc_max(i) = max(binall(disc_idx,2));
end

disc_all = [binmid',disc_mean',disc_min',disc_max'];
end

