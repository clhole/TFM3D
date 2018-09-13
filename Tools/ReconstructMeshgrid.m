function [ x, y, z ] = ReconstructMeshgrid( vect2D )
%RECONSTRUCTMESHGRID reconstructs a regular 3D grid in 'meshgrid' format
%                    from a 2D vector containing the respective coordinates
%Input:
%   <vect2D>         'meshgrid' coordinates in 2D vector format [x,y,z]
%Output:
%   <x>              x coordinates in 3D 'meshgrid' format
%   <y>              y coordinates in 3D 'meshgrid' format
%   <z>              z coordinates in 3D 'meshgrid' format
%CL

% Get directional distances between subsequent nodes
d = diff(vect2D);

% Get grid spacing and coordinate range
spacing = [max(d(:,1)),max(d(:,2)),max(d(:,3))];
range = [min(vect2D);max(vect2D)]';

% Determine size of the 3D 'meshgrid' matrix
matsize = round((range(:,2)-range(:,1))./spacing' + 1);

% Create 3D 'meshgrid' matrix
% x = reshape(vect2D(:,1),[matsize(2),matsize(1),matsize(3)]);
% y = reshape(vect2D(:,2),[matsize(2),matsize(1),matsize(3)]);
% z = reshape(vect2D(:,3),[matsize(2),matsize(1),matsize(3)]);

[x, y, z] = meshgrid(range(1,1):spacing(1):range(1,2),range(2,1):spacing(2):range(2,2),range(3,1):spacing(3):range(3,2));

end

