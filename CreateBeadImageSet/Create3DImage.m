function [ im ] = Create3DImage( vxsize, range, r, list, varargin )
%CREATE3DIMAGE creates a 3D bead image. Based on code by Tobias Schoch.
%              The following functions are used:
%              - gibsonlanni and shiftPSFAxial by Claude Holenstein
%              - convolution3D_FFTdomain by Christopher Coello
%Input:
%  <vxsize>    voxel size [dx,dy,dz]
%  <range>     coordinate range of the image volume [x1,x2;y1,y2;z1,z2]
%  <r>         bead radius
%  <list>      bead positions [x,y,z]
%  <PSFparam>  (optional) parameters for PSF calculation: [NA,ni,ns]
%Output:
%   <im>       3D bead image
%CL

% Define parameters
if ~isempty(varargin)
    usePSF = true; % PSF on
    PSFparam = varargin{1};
    %--------------------------------------------------------
    zp = 0;
    NA = PSFparam(1);
    ni = PSFparam(2);
    ns = PSFparam(3);
    wl = 610;
    microscope = 'Confocal';
    libpath = strcat(pwd,'\CreateBeadImageSet\PSF_Library');
    %--------------------------------------------------------
else
    usePSF = false; % PSF off
end

% Blue beads
color = 1;

% Prepare array for image
n = ceil((range(:,2)-range(:,1))./vxsize')+1;
im = zeros(n(1),n(2),n(3));

% Prepare scale
nbeads = size(list,1);

% Convert bead positions to pixel
list(:,1) = (list(:,1)-range(1,1))./vxsize(1);
list(:,2) = (list(:,2)-range(2,1))./vxsize(2);
list(:,3) = (list(:,3)-range(3,1))./vxsize(3);

% Calculate integer part
nlist = floor(list);

% Calculate half widths
hw = r*ones(3,1)./vxsize';
    
% Create sphere
nr = (round(hw./(.6))+2)*2+1;
m = zeros(nr(1),nr(2),nr(3));
m2 = zeros(nr(1),nr(2),nr(3));
c = (nr+1)./2;

if usePSF % Normal sphere for use with PSF
    in = 1;
    out = 1;
else
    in = 1-sqrt(3)./2;
    out = 1+sqrt(3)./2;
end

for i=1:nr(1)
    for j=1:nr(2)
        for k=1:nr(3)
            rpx = sqrt(((i-c(1))./hw(1)).^2+((j-c(2))./hw(2)).^2+((k-c(3))./hw(3)).^2);
            if rpx<=out
                if rpx<=in
                    m(i,j,k) = color;
                    m2(i,j,k) = color;
                else 
                    m(i,j,k) = (-1/sqrt(3)*rpx+(2+sqrt(3))/(2*sqrt(3)))*color;
                    m2(i,j,k) = (-1/sqrt(3)*rpx+(2+sqrt(3))/(2*sqrt(3)))*color;
                end
                if(any(hw<0.5))
                    m2(i,j,k)=m2(i,j,k)*prod(min(hw./0.5,[1;1;1]));
                end
            end
        end
    end
end

% Get PSF
if usePSF
    
    % Load psf library
    psfindex_file = load(strcat(libpath,'\psfindex.mat'));
    psfindex = psfindex_file.psfindex;
    npsf = size(psfindex,1);
    %--------------------------------------------------
    nzmax = 2*(range(3,2)-range(3,1))./vxsize(3);
    nzmax = 2*round((nzmax-1)/2)+1;
    np = max(nr.*5,[nr(1),nr(2),nzmax]');
    np(3) = min([np(3),2*round((10/vxsize(3)-1)/2)+1]);
    %--------------------------------------------------
    % Get psf parameter set
    psfparamset = [np',vxsize,zp,NA,ni,ns,wl];
    % Look for parameter set in library
    [~,psfnumber] = ismember(psfparamset,psfindex{:,1:11},'rows');
    if psfnumber ~= 0 && strcmp(psfindex.microscope{psfnumber},microscope)
        % Load existing psf from library
        fprintf(['Loading PSF %d from library... np=[%d,%d,%d], vxsize=',...
            '[%.3f,%.3f,%.3f], zp=%.2f, NA=%.2f, ni=%.2f, ns=%.2f, wl=%d, %s\n'],...
            psfnumber,psfindex{psfnumber,1:11},psfindex.microscope{psfnumber})
        psfdata = load(sprintf('%s\\psf%03d.mat',libpath,psfnumber));
        psf = psfdata.psf;
        np = psfdata.np;
    else % Calculate psf
        %------------------------------------------------------------
        psf = gibsonlanni(np',vxsize'*1e3,zp,NA,ni,ns,wl,microscope);
        %------------------------------------------------------------
        % Save to library
        psfindex{npsf+1,1:11} = psfparamset;
        psfindex.microscope{npsf+1} = microscope;
        save(strcat(libpath,'\psfindex.mat'),'psfindex');
        psfdata.psf = psf;
        psfdata.np = np;
        save(sprintf('%s\\psf%03d.mat',libpath,npsf+1),'-struct','psfdata');
    end
    psf_shift = shiftPSFAxial(psf,np);
    
    m2 = convolution3D_FFTdomain(psf_shift,m);
    m = m2./max(m2(:));
    
    nr = size(m)';
    c = (nr+1)./2;
    
end

% Shift positions to edge of window size
list = list-transpose((c-2)*ones(1,size(list,1)));
nlist = nlist-transpose((c-2)*ones(1,size(list,1)));

% Takes beads from list
for l = 1:nbeads
    
    % Take integer part of pos
    pos = nlist(l,:);
    
    % Take fractional part of pos
    shift = (list(l,:)-nlist(l,:));
    
    % Shift bead
    m_shift = ShiftImage(m,shift,[1,1,1]);
    
    left = max(ones(1,3),2-pos);
    right = min(nr,n-transpose(pos)+1);
    
%     for i = left(1):right(1)
%         for j = left(2):right(2)
%             for k = left(3):right(3)
%                 im(pos(1)+i-1,pos(2)+j-1,pos(3)+k-1) = im(pos(1)+i-1,pos(2)+j-1,pos(3)+k-1) + m_shift(i,j,k);
%             end
%         end
%     end
  
    im(pos(1)+left-1:pos(1)+right-1,pos(2)+left(2)-1:pos(2)+right(2)-1,pos(3)+left(3)-1:pos(3)+right(3)-1) =...
        im(pos(1)+left-1:pos(1)+right-1,pos(2)+left(2)-1:pos(2)+right(2)-1,pos(3)+left(3)-1:pos(3)+right(3)-1) + ...
        + m_shift(left(1):right(1),left(2):right(2),left(3):right(3));
    
end

end

% ShiftImage
%**************************************************************************
function [ m_shift ] = ShiftImage( m, u, pxdelta )
%SHIFTIMAGE shifts a 3D images <m> by the vector <u> in units of <pxdelta>
%           (<u> [um] <--> <pxdelta> [um/px] in each dimension)

% Convert to pixel
u = u./pxdelta;

% Take integer part of pos
pos = floor(u);

% Do integer shifting
m = circshift(m,pos);

% Take fractional part of pos
shift = u-pos;

% Shift sphere about the fractional part
% Prepare filters
firX = [0;1-shift(1);shift(1)];
firY = [0,1-shift(2),shift(2)];
firZ(:,:,1) = 0;
firZ(:,:,2) = 1 - shift(3);
firZ(:,:,3) = shift(3);

% Shift bead
m_shift = convn(m,firX,'same');
m_shift = convn(m_shift,firY,'same');
m_shift = convn(m_shift,firZ,'same');

end
