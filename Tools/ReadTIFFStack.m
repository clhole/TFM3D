function [ image ] = ReadTIFFStack( filename )
%READTIFFSTACK   reads a .tif image stack in 3D matrix format
%Input:
%   <filename>   image stack (.tif)
%Output:
%   <image>      3D matrix
%CL

% Get file information
imginfo = imfinfo(filename);
m = imginfo(1).Width;
n = imginfo(1).Height;
nimg = length(imginfo);

image = zeros(n,m,nimg,'uint16');

% Read data
tiffobj = Tiff(filename,'r');
for i = 1:nimg
   tiffobj.setDirectory(i);
   image(:,:,i) = tiffobj.read();
end
tiffobj.close();

% Convert to double
image = im2double(image);

end

