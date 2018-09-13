function stl2fv2sat( varargin )
%STL2FV2SAT      converts an .stl geometry to .mat patch and .sat format.
%                The geometry can be shifted to a desired location.
%                The following functions by Adam H. Aitkenhead are used:
%                - READ_stl
%                - WRITE_sat
%                - PLOT_3D_stl_patch
%Input:
%   <filename>   (optional) geometry file in .stl format
%   <loc>        (optional) new location of the vertex centroid [x,y,z]
%CL

% Parse input
switch nargin
    case 2
        [path,name,~] = fileparts(varargin{1});
        loc = varargin{2};
        gui = false;
    case 1
        [path,name,~] = fileparts(varargin{1});
        loc = [];
        gui = false;
    otherwise
        [filename,path] = uigetfile('*.stl');
        [~,name,~] = fileparts(filename);
        gui = true;
end

% Read .stl file
stlfile = fullfile(path,strcat(name,'.stl'));
[fvmat,normals] = READ_stl(stlfile,'auto');

% Remove any facets containing a zero-length edge
zeroarealist = max([min(fvmat(:,:,1)==fvmat(:,:,2),[],2),...
                    min(fvmat(:,:,2)==fvmat(:,:,3),[],2),...
                    min(fvmat(:,:,1)==fvmat(:,:,3),[],2)],[],2);
fvmat = fvmat(zeroarealist==0,:,:);
normals  = normals(zeroarealist==0,:,:);

% Calculate face centroids
facecent = (fvmat(:,:,1)+fvmat(:,:,2)+fvmat(:,:,3))/3;
% Calculate face areas
fareas = 1/2*sqrt(sum(cross(fvmat(:,:,2)-fvmat(:,:,1),fvmat(:,:,3)-fvmat(:,:,1)).^2,2));
% Calculate area weighted centroid
centroid(:,1) = sum(facecent(:,1).*fareas)/sum(fareas);
centroid(:,2) = sum(facecent(:,2).*fareas)/sum(fareas);
centroid(:,3) = sum(facecent(:,3).*fareas)/sum(fareas);

% Shift coordinates
if gui
    % Get new centroid coordinates
    answer = inputdlg({'x','y','z'},'Centroid',[1,40;1,40;1,40],...
        {num2str(centroid(1),9),num2str(centroid(2),9),num2str(centroid(3),9)});
    loc = [str2double(answer{1}),str2double(answer{2}),str2double(answer{3})];
     % Avoid shifting error due to num2str precision
    if max(abs(centroid-loc)) < 10^-6, loc = []; end;
end
if ~isempty(loc)
   nfaces = size(fvmat,1);
   fvmat(:,:,1) = fvmat(:,:,1) - ones(nfaces,1)*(centroid + loc);
   fvmat(:,:,2) = fvmat(:,:,2) - ones(nfaces,1)*(centroid + loc);
   fvmat(:,:,3) = fvmat(:,:,3) - ones(nfaces,1)*(centroid + loc); 
   facecent = facecent - ones(nfaces,1)*(centroid + loc);  
end

% Get all vertices for .mat polygon format
v_all = [fvmat(:,:,1);fvmat(:,:,2);fvmat(:,:,3)];
% Remove duplicate vertices
v = unique(v_all,'rows');
% Get new vertex indices for each face
nfaces = size(fvmat,1);
f = zeros(nfaces,3);
for n = 1:nfaces
   [~,f(n,1)] = ismember(fvmat(n,:,1),v,'rows');
   [~,f(n,2)] = ismember(fvmat(n,:,2),v,'rows');
   [~,f(n,3)] = ismember(fvmat(n,:,3),v,'rows');
end 
UserData.facecent = facecent;
UserData.normals = normals;
fvcell = struct('faces',f,'vertices',v,'UserData',UserData);
% Write .mat file
matfile = fullfile(path,strcat(name,'.mat'));
save(matfile,'-struct','fvcell');

% Write .sat file
satfile = fullfile(path,strcat(name,'.sat'));
errmsg = WRITE_sat(satfile,fvmat,normals);
if sum(errmsg) > 0
  fprintf(' File conversion failed.\n')
else
  fprintf(' File saved as %s\n',satfile)
end

% Plot geometry
figure
patch(fvcell,'FaceColor','g','FaceAlpha',0.2)
xlabel('x [µm]')
ylabel('y [µm]')
zlabel('z [µm]')
axis equal
grid on
view(3)

end