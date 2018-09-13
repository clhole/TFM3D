function [ nodes, u ] = GDVC( img0, img1, rbead, vxsize, range, maxwinsize, minwinsize, spacing )
%GDVC           Grid-based DVC to calculate displacements between two 3D images
%Input:
%  <img0>       bead image at t=t0 (with cell)
%  <img1>       bead image at t=t1 (without cell) 
%  <rbead>      bead radius [µm]
%  <vxsize>     voxel size [dx;dy;dz] [µm]
%  <range>      nodal range [x1,x2;y1,y2;z1,z2] [µm]
%  <maxwinsize>	initial search window size [px]
%  <minwinsize>	minimum search window size [px]
%  <spacing>    search window spacing [px]
%Output:
%  <nodes>      grid nodes where displacements were calculated [x,y,z] [µm]
%  <u>          grid node displacements [ux,uy,uz] [µm]
%CL

% Get parameters [px]
imgsize = size(img0);
rbead_px = rbead./vxsize;

% Prepare grid
for d = 1:3
    gv{d} = 1:spacing(d):imgsize(d);
end
[x,y,z] = meshgrid(gv{1},gv{2},gv{3});
nodes = [x(:),y(:),z(:)];
nnodes = size(nodes,1);

% Initialization
u = zeros(nnodes,3);
nerrors = 0;
ind_del = [];
clearmsg = '';

% Go through all grid nodes
for i = 1:nnodes
    
    % Set starting search window size
    winsize = maxwinsize;

    % Reduce window size
    while all(winsize>minwinsize)
        
        % Get window of the undeformed state image
        window1 = GetWindowAt(img0,winsize,nodes(i,:));
        % Get search window of the deformed state image
        window2 = GetWindowAt(img1,winsize,nodes(i,:)+round(u(i,:)));

        % Correllation fit
        [u_win,exitflags] = GaussFit(window1,window2,rbead_px);
        u(i,:) = round(u(i,:)) + u_win;
        
        %Check for bad fits or impossible displacements
        if exitflags~=3 || any(abs(u_win)>=floor(winsize/2)) || any(nodes(i,:)+round(u(i,:))<-1)
            nerrors = nerrors + 1;
            ind_del(nerrors) = i;
            break;
        end
        
        % Reduce the search window size
        winsize = round(winsize*2/3);

    end
    
    msg = sprintf('Progress: %.2f%%%%. Number of errors: %d',i/nnodes*100,nerrors);
    fprintf([clearmsg,msg])
    clearmsg = repmat(sprintf('\b'),1,length(msg)-1);
    
end

% Remove bad fits
nodes(ind_del,:) = [];
u(ind_del,:) = [];

% Convert nodal coordinates and displacements to [µm]
nodes(:,1) = (nodes(:,1)-1)*vxsize(1)+range(1,1);
nodes(:,2) = (nodes(:,2)-1)*vxsize(2)+range(2,1);
nodes(:,3) = (nodes(:,3)-1)*vxsize(3)+range(3,1);
u(:,1) = u(:,1)*vxsize(1);
u(:,2) = u(:,2)*vxsize(2);
u(:,3) = u(:,3)*vxsize(3);

end