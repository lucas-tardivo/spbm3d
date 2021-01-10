function [x] = addBorders(i, bordersize)

[nrows,ncols] = size(i);
x = zeros(nrows+(bordersize*2),ncols+(bordersize*2));
x(:,:) = 0.0;
x((bordersize+1):end-bordersize,(bordersize+1):end-bordersize,:) = i;

end