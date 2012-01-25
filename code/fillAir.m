function [ output ] = fillAir( tissue )
% fillAir() fills gaps in data with air
%   Once you import all of the data using loadNII(), run it thought this to
%   fill in the remaining spaces with air.  

airdata = [1 0 0 1006 1.3 0.026 0];

%  Picks out air spots
a = find(tissue(:,:,:,1) == 0);  
[x y z t] = ind2sub(size(tissue),a);

for i = 1:length(a)
    tissue(x(i),y(i),z(i),:) = airdata;
end

output = tissue;

end