function InpolyhedronTest()
%INPOLYHEDRONTEST is a testing tool for the function <inpolyhedron.m>
%CL

% Parameters
%**************************************************************************
filename_cell = 'Neuron.msh';       % Cell model
nnodes = [100,100,1];               % Node numbers [nx,ny,nz]
%**************************************************************************

% Load cell model
fv = ImportMesh(filename_cell,true);
% Change face normals form inward to outward pointing
fv.faces = fliplr(fv.faces);

% Create grid
minv = min(fv.vertices);
maxv = max(fv.vertices);
gridx = minv(1):(maxv(1)-minv(1))/(nnodes(1)-1):maxv(1);
gridy = minv(2):(maxv(2)-minv(2))/(nnodes(2)-1):maxv(2);
if nnodes(3) ~= 1
    gridz = minv(3):(maxv(3)-minv(3)/(nnodes(3)-1)):maxv(3);
else
    gridz = (maxv(3)+minv(3))/2;
end
[x,y,z] = meshgrid(gridx,gridy,gridz);

% Check if grid points are inside the object
in = inpolyhedron(fv,gridx,gridy,gridz);

% Calculate face centers and face normals
normals = zeros(size(fv.faces,1),3);
centers = zeros(size(fv.faces,1),3);
for i = 1:size(fv.faces,1)
    normals(i,:) = cross((fv.vertices(fv.faces(i,2),:)-fv.vertices(fv.faces(i,1),:)),...
        (fv.vertices(fv.faces(i,3),:)-fv.vertices(fv.faces(i,2),:)));
    centers(i,1) = (fv.vertices(fv.faces(i,1),1)+fv.vertices(fv.faces(i,2),1)+...
        fv.vertices(fv.faces(i,3),1))/3;
    centers(i,2) = (fv.vertices(fv.faces(i,1),2)+fv.vertices(fv.faces(i,2),2)+...
        fv.vertices(fv.faces(i,3),2))/3;
    centers(i,3) = (fv.vertices(fv.faces(i,1),3)+fv.vertices(fv.faces(i,2),3)+...
        fv.vertices(fv.faces(i,3),3))/3;
end

figure
hold on
% Plot cell and face normals
patch(fv,'FaceColor','g','FaceAlpha',0.4)
quiver3(centers(:,1),centers(:,2),centers(:,3),normals(:,1),normals(:,2),normals(:,3))
% Plot points inside (blue) and outside (red) the cell
plot3(x(in), y(in), z(in),'b+')
plot3(x(~in),y(~in),z(~in),'r+')
% Set figure properties
xlabel('x [µm]')
ylabel('y [µm]')
zlabel('z [µm]')
axis equal
view(3)

end

