function WriteTIFFStack( mat )
%WRITETIFFSTACK  saves a 3D matrix as a 16 bit .tif file
%Input:
%   <mat>        3D matrix
%CL

% Conversion to 16 bit format
mat = im2uint16(mat);

% Write data
imwrite(mat(:,:,1),'stack.tif')
for k = 2:size(mat,3)
    imwrite(mat(:,:,k),'stack.tif','writemode','append');
end

end