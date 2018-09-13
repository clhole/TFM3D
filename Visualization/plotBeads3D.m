% Tabula Rasa
clear all
close all

% Image generation parameters for 3D TFM simulation (for mode = 'sim*')
%---------------------------------------------------------  ----------------------
folder = '3_Sphere20um_x3';         % Project folder in 'TFM 3D\FEAData'
cbeads = [0.3];                   % Bead concentration [1/µm^3]       [1xn]
rbead = [0.1];                     % Bead radius [µm]              [1xn]
vxsize = [0.1,0.1,0.1];            % Voxel size [µm]                   [nx3]
PSFparam = [1.2,1.33,1.33];        % PSF: [NA,ni,ns], [] = PSF off     [nx3]

% Define ROI
% ROI = ones(3,1)*[-15 15];
ROI = [-15 15;-0.2 0.2;-15 15];

% dataset = CreateBeadImageSet(folder,cbeads(i),rbead,vxsize,PSFparam);

% Import cell surface patch from .mat file
file_cell = strcat(pwd,'\FEAData\',folder,'\',folder,'.mat');
fvcell = load(file_cell);

% Import FEA results from ANSYS Workbench (as interpolants)
file_wb = strcat(pwd,'\FEAData\',folder,'\',folder,'.wbpj');
[fuFEA,fTFEA,frFEA,range] = ImportFEAResults(file_wb);

% Generate bead positions at t=t0 (with the cell)
beadpos0 = RandBeadPos(range,cbeads,max(rbead),fvcell);
nbeads = size(beadpos0,1);

% Determine displacements and displaced bead positions at t=t1 (without the cell)
uFEA = zeros(nbeads,3);
beadpos1 = zeros(nbeads,3);
dim = ['x','y','z'];
for d = 1:3
    uFEA(:,d) = fuFEA.(dim(d))(beadpos0);
    beadpos1(:,d) = beadpos0(:,d) + uFEA(:,d);
end

% Create image pair
PSFparamset = PSFparam();
img0 = Create3DImage(vxsize,range,rbead,beadpos0,PSFparamset);
img1 = Create3DImage(vxsize,range,rbead,beadpos1,PSFparamset);

vertices = fvcell.vertices;
faces = fvcell.faces;

% Get beads within the region of interest
ROI_px(:,1) = ceil(bsxfun(@rdivide,ROI(:,1)-range(:,1),vxsize')); % [px]
ROI_px(:,2) = ceil(bsxfun(@rdivide,ROI(:,2)-range(:,1),vxsize')); % [px]
for d = 1:3;
    subvol{d} = ROI_px(d,1)+1:ROI_px(d,2)+1;
end;
img0 = img0(subvol{1},subvol{2},subvol{3});
img1 = img1(subvol{1},subvol{2},subvol{3});
beadpos0 = beadpos0(beadpos0(:,1)>=ROI(1,1) & beadpos0(:,1)<=ROI(1,2) &...
                    beadpos0(:,2)>=ROI(2,1) & beadpos0(:,2)<=ROI(2,2) &...
                    beadpos0(:,3)>=ROI(3,1) & beadpos0(:,3)<=ROI(3,2),:);
imgsize = size(img0);

% Convert bead positions to [px] (center pixel coordinates)
beadpos0_px(:,1) = (beadpos0(:,1)-range(1,1))./vxsize(1); % [px]
beadpos0_px(:,2) = (beadpos0(:,2)-range(2,1))./vxsize(2); % [px]
beadpos0_px(:,3) = (beadpos0(:,3)-range(3,1))./vxsize(3); % [px]

% Convert cell model to [px] (center pixel coordinates)
fvcell_px.vertices(:,1) = (vertices(:,1)-range(1,1))./vxsize(1); % [px]
fvcell_px.vertices(:,2) = (vertices(:,2)-range(2,1))./vxsize(2); % [px]
fvcell_px.vertices(:,3) = (vertices(:,3)-range(3,1))./vxsize(3); % [px]
fvcell_px.faces = faces;

% Prepare grid in [px]
[x,y,z] = ndgrid(ROI_px(1,1):ROI_px(1,2),ROI_px(2,1):ROI_px(2,2),ROI_px(3,1):ROI_px(3,2));
% Intensity threshold to display beads
tsh = max(img0(:))/8;

figure(1)
hold on
patch(isosurface(x,y,z,img0,tsh),'FaceColor','r','FaceAlpha',1,'EdgeAlpha',0.1);
patch(isosurface(x,y,z,img1,tsh),'FaceColor','g','FaceAlpha',1,'EdgeAlpha',0.1);
daspect([1,1,1])
view([0 0]); axis tight
camlight 
lighting gouraud
% Plot cell
patch(fvcell_px,'FaceColor','g','FaceAlpha',1)
set(gca,'XTick', [250 400 550],'XTickLabel',{'-15','0','15'},'FontSize',20)
set(gca,'YTick', [250 400 550],'YTickLabel',{'-15','0','15'},'FontSize',20)
set(gca,'ZTick', [250 400 550],'ZTickLabel',{'-15','0','15'},'FontSize',20)
box on, grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
print(strcat(pwd,'\figures\cell_beads_3D'),'-dtiff','-r0')

% Second image with beads as bitmap (on/off)
sx = size(img0,1);
sy = size(img0,2);
sz = size(img0,3);

img_comp= uint8(zeros(sx,sy,3));  
img_comp(:,:,1) = im2uint8(img0(:,:,150));
img_comp(:,:,2) = im2uint8(img1(:,:,150));

img_comp2= uint8(zeros(sx,sz,3));  
img_comp2(:,:,1) = im2uint8(squeeze(img0(:,ceil(sy/2),:)));
img_comp2(:,:,2) = im2uint8(squeeze(img1(:,ceil(sy/2),:)));
comp_fin = permute(img_comp2,[2 1 3]);

figure(2)
% Plot cell
patch(fvcell_px,'FaceColor','none','FaceAlpha',1,'EdgeColor','white')
hold on
% th = 0:pi/50:2*pi;
% r = 10./vxsize(1);
% x0 = 400;
% z0 = x0;
% xunit = r * cos(th) + x0;
% zunit = r * sin(th) + z0;
% h = plot(xunit, zunit);

xImage = [250 550; 250 550];   %# The x data for the image corners
yImage = [400 400; 400 400];   %# The y data for the image corners
zImage = [250 250; 550 550];   %# The z data for the image corners
surf(xImage,yImage,zImage,...    %# Plot the surface
     'CData',comp_fin,...
     'FaceColor','texturemap');
view([0 0]); axis tight square
% camlight 
% lighting gouraud
set(gca,'XTick', [250 400 550],'XTickLabel',{'-15','0','15'},'FontSize',20)
set(gca,'YTick', [250 400 550],'YTickLabel',{'-15','0','15'},'FontSize',20)
set(gca,'ZTick', [250 400 550],'ZTickLabel',{'-15','0','15'},'FontSize',20)
box on, grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off

print(strcat(pwd,'\figures\beads_bitmap'),'-dtiff','-r0')

[x2, y2, z2] = meshgrid(ROI_px(1,1):ROI_px(1,2),ROI_px(2,1):ROI_px(2,2),ROI_px(3,1):ROI_px(3,2));

figure(3)
h = slice(x2,y2,z2,permute(img0, [2 1 3]),[],[400],[]);
set(h, 'EdgeColor','none', 'FaceColor','interp')
alpha(.8)
colormap(jet)
daspect([1,1,1])
view([0 0]); axis tight


h(1).EdgeColor = 'none';


% h(2).EdgeColor = 'none';
% function h = circle(x,y,r)
%     hold on
%     th = 0:pi/50:2*pi;
%     xunit = r * cos(th) + x;
%     yunit = r * sin(th) + y;
%     h = plot(xunit, yunit);
%     hold off
% end

