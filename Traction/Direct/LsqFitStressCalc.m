function [ spos, s ] = LsqFitStressCalc( upos, ucalc, E, nu, spacing, winsize )
%LSQFITSTRESSCALC directly calculates local LSQ fitted displacement field
%                 gradients to obtain material strains which are then used
%                 to compute stresses. Based on code by Tobias Schoch.
%Input:
%	<upos>        positions where displacements are known [µm]
%	<ucalc>       known displacements [µm]
%	<E>           Young's modulus of the substrate [Pa]
%	<nu>          Poisson's ratio of the substrate
%   <spacing>     grid spacing for grad(u) calculation [µm]
%   <winsize>     window size for grad(u) averaging [grid units]
%Output:
%   <spos>        positions where stress tensors are calculated [µm]
%	<s>           stress tensors in Nye notation [sx,sy,sz,sxy,sxz,syz] [Pa]
%CL

%**************************************************************************
% Define smoothing regulation parameters
p = [.02 .02 .01];
a = (log(p(1:2))-log(p(3)))./log(0.1);
%**************************************************************************

% Determine Lamé constants (homogeneous isotropic linear elastic material)
lambda = nu*E/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));

% Slicewise 2D spline interpolation of the displacement field
%**************************************************************************

% Prepare 2D grid within the range of displacement data
xran = min(upos(:,1)):spacing:max(upos(:,1));
yran = min(upos(:,2)):spacing:max(upos(:,2));
nx = size(xran,2);
ny = size(yran,2);
[x,y] = ndgrid(xran,yran);

% Get intrinsic slices
zran = unique(upos(:,3));
nslices = size(zran,1);

u_itp = zeros(nx,ny,nslices,3);
clearmsg = '';

for n = 1:nslices
    
    % Find points on slice where displacements are known
    ind = find(upos(:,3)==zran(n));
    
    % Calculate oversampling factor and regulate smoothness
    f = length(ind)./(nx*ny);
    p_reg = p(1:2)./(f.^a);
    p_reg = min(p_reg,[1,1]);
    p_reg(3) = p(3);
    
    for d = 1:3
        
        % Define fitting options
        ft = fittype('loess');
        opts = fitoptions(ft);
        opts.Robust = 'off';
        opts.Span = p_reg(d);
        opts.Normalize = 'on';
        
        % Fitting
        fit_obj = fit(upos(ind,1:2),ucalc(ind,d),ft,opts);        
        u_itp(:,:,n,d) = fit_obj(x,y);
        
    end

    msg = sprintf('Slicewise displacement field interpolation... %d/%d',n,nslices);
    fprintf('%s',[clearmsg,msg])
    clearmsg = repmat(sprintf('\b'),1,length(msg));
    
end

% Smoothing of the interpolated displacement field in z dir. using splines
%**************************************************************************

for indx = 1:nx
    for indy = 1:ny
        for d = 1:3
        
            % Get z data
            zdata = u_itp(indx,indy,:,d);
            zran_temp = zran;

            % Remove NaNs
            ind = find(isnan(zdata));
            zdata(ind) = [];
            zran_temp(ind) = [];
            if size(zdata,1) < 10
                continue;
            end

            % Fitting
            fit_obj = csaps(zran_temp',zdata,p);
            u_itp(indx,indy,:,d) = fnval(fit_obj,zran_temp);

        end
    end
end

% Strain and stress calculation
%**************************************************************************

s = zeros(nx*ny*nslices,6);
spos = zeros(nx*ny*nslices,3);
n = 0;
clearmsg = '';
fprintf('\n')

for indx = 1:nx
    for indy = 1:ny
        for indz = 1:nslices
            
            n = n+1;
            
            % Gradient calculated over a (nxnxn)-window by a linear lsq fit
            Du = GradientAt(u_itp,[spacing;spacing;spacing],winsize,[indx,indy,indz]);
            
            % Calculate strain tensor
            epsilon = 1/2*(Du + Du');
            
            % Calculate stress tensor
            sigma = lambda*trace(epsilon)*eye(3) + 2*mu*epsilon;
            
            % Store resulting state of stress in vector format: [sx,sy,sz,sxy,sxz,syz]
            s(n,:) = [sigma(1,1),sigma(2,2),sigma(3,3),sigma(1,2),sigma(1,3),sigma(2,3)];
            
            % Store nodal coordinates
            spos(n,:) = [xran(indx),yran(indy),zran(indz)];
            
        end
    end
    
    msg = sprintf('Stress field calculation... %.0f%%',indx/nx*100);
    fprintf('%s',[clearmsg,msg])
    clearmsg = repmat(sprintf('\b'),1,length(msg));
    
end

fprintf('\n')

end