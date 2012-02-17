function [ array_out ] = lambw( array_in )
% lambw Wrapper for wapr()
% Available:  http://www.mathworks.com/matlabcentral/fileexchange/3644-real-values-of-the-lambert-w-function/content/Lambert/wapr.m
%   Dwapr() doesn't work any arrays over Nx1, so this steps through the
%   full matrix and gives the rows to wapr.  Works pretty fast.
%#codegen

if ndims(array_in) ~= 3
    error('This only works (for now) with a three dimensional array.')
end

xmax = size(array_in,1);
ymax = size(array_in,2);

array_out = zeros(size(array_in));
  for ix=1:xmax
      for iy=1:ymax
          array_out(ix,iy,:) = wapr(array_in(ix,iy,:));
      end
  end
end