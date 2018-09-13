function [ window, img_padded, ind_win, padsize ] = GetWindowAt( img, winsize, ind )
%GETWINDOWAT    extracts a window of a 2D or 3D image, centered around the 
%               pixel/voxel <ind>
%Input
%   <img>       2D or 3D array
%   <winsize>   window size [lx,ly,(lz)]
%   <ind>       index of the center voxel [x,y,(z)]
%Output
%   <window>    2D or 3D window
%CL

% Padding
padsize = floor(winsize/2)+2;
img_padded = padarray(img,padsize);

% Extract the window
winrange = [ind+2;ind+winsize+1]';
for d = 1:ndims(img)
    ind_win{d} = winrange(d,1):winrange(d,2);
end
window = img_padded(ind_win{:});

end