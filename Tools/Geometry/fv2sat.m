function fv2sat( fvpatch, filename_sat)
%FV2SAT             converts a patch class object to .sat format using the
%                   function WRITE_sat by Adam H. Aitkenhead.
%Input:
%  <fvpatch>        patch class object
%  <filename_sat>   filename of the .sat file
%CL

% Get faces and vertices
faces = fvpatch.faces;
vertices = fvpatch.vertices;

% Create vertex matrix for WRITE_sat
nfaces = size(faces,1);
fvmat = zeros(nfaces,3,3);
for i = 1:nfaces
    fvmat(i,:,1) = vertices(faces(i,1),:);
    fvmat(i,:,2) = vertices(faces(i,2),:);
    fvmat(i,:,3) = vertices(faces(i,3),:);
end

% Calculate face normals
normals = zeros(nfaces,3);
for i = 1:nfaces
    crossp = cross((vertices(faces(i,2),:)-vertices(faces(i,1),:)),...
        (vertices(faces(i,3),:)-vertices(faces(i,2),:)));
    normals(i,:) = crossp/norm(crossp);
end

% Write .sat file
errmsg = WRITE_sat(filename_sat,fvmat,normals);
if sum(errmsg) > 0
  fprintf(' File conversion failed.\n')
else
  fprintf(' File saved as %s\n',filename_sat)
end

end