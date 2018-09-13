function PlotBeadImagePair( imgpair, varargin )
%PLOTBEADIMAGEPAIR   displays the generated 3D bead image pair
%Input:
%   <imgpair>        structure array containing the fields
%                    .vxsize       voxel size [µm]
%                    .range        nodal range [x1,x2;y1,y2;z1,z2] [µm]
%                    .fvcell       cell model patch object [µm]
%                    .beadpos0     given neutral bead positions [µm]
%                    .img0         undeformed state bead image
%                    .img1         deformed state bead image
%   <ROI>            (optional) region of interest [µm]
%   <datasetnumber>  (optional) datasetnumber to display in the figure title
%CL

% Get images and parameters
vxsize = imgpair.imgparam.vxsize;
range = imgpair.imgparam.range;
vertices = imgpair.imgparam.fvcell.vertices;
faces = imgpair.imgparam.fvcell.faces;
beadpos0 = imgpair.imgparam.beadpos0;

% Get optional input
if nargin == 3
    if size(varargin{1},2)>1
        ROI = varargin{1};
        datasetnumber = varargin{2};
    else
        ROI = varargin{2};
        datasetnumber = varargin{1};
    end
elseif nargin == 2
    if size(varargin{1},2)>1
        ROI = varargin{1};
        datasetnumber = 0;
    else
        ROI = range;  
        datasetnumber = varargin{1};
    end
else
    ROI = range;
    datasetnumber = 0;
end
% Interpolation on cell vertices accepts no ROI
if strcmp('cell vertices',ROI)
    ROI = range;
end
% Display the number of the dataset in the figure title
if datasetnumber ~=0
    reprdataset = sprintf('Representative dataset: No. %03d',datasetnumber);
else
    reprdataset = '';
end

% Get beads within the region of interest
ROI_px(:,1) = ceil(bsxfun(@rdivide,ROI(:,1)-range(:,1),vxsize')); % [px]
ROI_px(:,2) = ceil(bsxfun(@rdivide,ROI(:,2)-range(:,1),vxsize')); % [px]
for d = 1:3;
    subvol{d} = ROI_px(d,1)+1:ROI_px(d,2)+1;
end;
img0 = imgpair.img0(subvol{1},subvol{2},subvol{3});
img1 = imgpair.img1(subvol{1},subvol{2},subvol{3});
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
tsh = max(img0(:))/4;

figure
hold on
% Plot beads
scatter3(beadpos0_px(:,1),beadpos0_px(:,2),beadpos0_px(:,3),'r+')
patch(isosurface(x,y,z,img0,tsh),'FaceColor','b','FaceAlpha',0.3);
patch(isosurface(x,y,z,img1,tsh),'FaceColor','r','FaceAlpha',0.3);
legend('Generated bead positions at t=t0','Generated beads at t=t0','Generated beads at t=t1')
% Plot cell
patch(fvcell_px,'FaceColor','g','FaceAlpha',0.3)
% Set figure properties
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.7 0.8])
title({sprintf('%s',reprdataset);'Generated bead distribution';strrep(sprintf(['%dx%dx%d voxels, '...
    '[%.1f,%.1f; %.1f,%.1f; %.1f,%.1f] µm, %d beads'],imgsize',ROI',size(beadpos0_px,1)),'-0','0')})
xlabel('x [voxel]')
ylabel('y [voxel]')
zlabel('z [voxel]')
axis equal
axis([ROI_px(1,:),ROI_px(2,:),ROI_px(3,:)])
grid on
view(3)

end

