function [ img_filt ] = LocThrFilter( img, filtwin )
%LOCTHRFILTER    filters background noise by local thresholding.
%                Inspired by the work of James H-C. Wang, 2006.
%Input:
%  <img>         2D or 3D image in matrix format
%  <filtwin>     size of the subwindows for local thresholding
%Output:
%  <image_filt>  filtered 2D or 3D image in matrix format
%CL

% Image size and dimensions
imgsize = size(img);
n = ndims(img);

% Divide the image into subwindows
for dim = 1:n
    % Get the number of subwindows of the given size and specify distances
    nsub_int(dim) = floor(imgsize(dim)/filtwin(dim));
    dist{dim} = filtwin(dim)*ones(1,nsub_int(dim));
    % Determine the margin size
    margin(dim) = rem(imgsize(dim),filtwin(dim));
    if margin(dim) ~= 0
        % Add distance of margin subwindow
        dist{dim}(1,end+1) = margin(dim);
    end
end
% Divide the image
if n == 3
    subwindows = mat2cell(img,dist{1},dist{2},dist{3});
else
    subwindows = mat2cell(img,dist{1},dist{2});
end
nsub = numel(subwindows);

% Filtering
for n = 1:nsub
    for rep = 1:4
        thr = mean(subwindows{n}(:))/100;
        subwindows{n} = subwindows{n}.^2;
        subwindows{n}(subwindows{n}<=thr) = 0;
    end
end

% Reassemble the image
img_filt = cell2mat(subwindows);

end

