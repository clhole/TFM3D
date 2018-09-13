function [ beadpos, u ] = PDVC( img0, img1, rbead, vxsize, range, umax, cbeads_win )
%BDVC            Particle-based DVC to calculate displacements between two 3D images.
%                Based on code by Tobias Schoch, 2012.
%Input:
%  <img0>        bead image at t=t0 (with cell)
%  <img1>        bead image at t=t1 (without cell) 
%  <rbead>       bead radius [µm]
%  <vxsize>      voxel size [dx;dy;dz] [µm]
%  <range>       nodal range [x1,x2;y1,y2;z1,z2] [µm]
%  <umax>        expected maximum displacement [µm]
%  <cbeads_win>  average number of beads per search window
%Output:
%  <beadpos>     detected bead positions [x,y,z] [µm] at t=t0
%  <u>           bead displacements in the format [ux,uy,uz] [µm]
%CL

% Get parameters [px]
imgsize = size(img0);
rbead_px = rbead./vxsize;
umax_px = umax./vxsize;

% Find beads (local maximum)
%**************************************************************************
fprintf('Locating beads... ')
beadpos = zeros(1,3);
nbeads = 0;

% Local threshold filtering using 7x7x7 voxel boxes
img0_filt = LocThrFilter(img0,[7,7,7]);

% Find beads
clearmsg = '';
for x = 2:imgsize(1)-1
    for y = 2:imgsize(2)-1 
        for z = 2:imgsize(3)-1
                              
            % Create a 3x3x3 voxel box around the current voxel
            box = img0_filt(x-1:x+1,y-1:y+1,z-1:z+1);
            % Check if there is a maximum intensity
            if ~all(box(:)==box(1))
                % Check if the maximum intensity is at the current voxel
                if box(2,2,2) == max(box(:))
                    nbeads = nbeads+1;
                    % Bead positions in [px], indexing starts at [1,1,1]
                    beadpos(nbeads,:) = [x,y,z];
                end
            end

        end
    end
    
    msg = sprintf('%.0f%% searched, %d beads found',x/(imgsize(1)-1)*100,nbeads);
    fprintf('%s',[clearmsg,msg])
    clearmsg = repmat(sprintf('\b'),1,length(msg));
    
end

% Bead tracking and bead displacement calculation (DVC method)
%**************************************************************************
fprintf('\nTracking beads... ')

% Determine maximum search window size depending on maximum expected displacement
maxwinsize = min(round(2.7*umax_px),imgsize);
% Determine minimum search window size from the maximum of minimum distance
% between beads and average distance between beads
dbeads_min = round(4*rbead_px);
dbeads_avg = round((prod(imgsize)/nbeads*cbeads_win)^(1/3))*ones(1,3);
minwinsize = max(dbeads_min,dbeads_avg);

% Initialization
u = zeros(nbeads,3);
nerrors = 0;
ind_del = [];
clearmsg = '';

% Go through all beads
for i=1:nbeads

    % Set starting search window size
    winsize = maxwinsize;
    
    % Reduce window size
    while all(winsize>minwinsize)
        
        % Get search window at undisplaced bead position
        window1 = GetWindowAt(img0,winsize,beadpos(i,:));
        % Get search window at estimated displaced bead position
        window2 = GetWindowAt(img1,winsize,beadpos(i,:)+round(u(i,:)));
        
        % Correllation fit
        [u_win,exitflags] = GaussFit(window1,window2,rbead_px);
        u(i,:) = round(u(i,:)) + u_win;

        %Check for bad fits or impossible displacements
        if exitflags~=3 || any(abs(u_win)>=floor(winsize/2)) || any(beadpos(i,:)+round(u(i,:))<-1)
            nerrors = nerrors + 1;
            ind_del(nerrors) = i;
            break;
        end
    
        % Reduce the search window size
        winsize = round(winsize*2/3);

    end
    
    msg = sprintf('%d. Number of errors: %d',i,nerrors);
    fprintf([clearmsg,msg])
    clearmsg = repmat(sprintf('\b'),1,length(msg));

end

fprintf('\n')

% Remove bad fits
beadpos(ind_del,:) = [];
u(ind_del,:) = [];

% Convert bead positions and displacements to [µm]
beadpos(:,1) = (beadpos(:,1)-1)*vxsize(1)+range(1,1);
beadpos(:,2) = (beadpos(:,2)-1)*vxsize(2)+range(2,1);
beadpos(:,3) = (beadpos(:,3)-1)*vxsize(3)+range(3,1);
u(:,1) = u(:,1)*vxsize(1);
u(:,2) = u(:,2)*vxsize(2);
u(:,3) = u(:,3)*vxsize(3);

end