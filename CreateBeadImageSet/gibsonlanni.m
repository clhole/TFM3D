function out = gibsonlanni(varargin)

fprintf('Calculating PSF...    ')

% Input parser
ip = inputParser;
addOptional(ip, 'n', [128,128,129]);
addOptional(ip, 'res', [50,50,100]);
addOptional(ip, 'zp', 0, @(x)isnumeric(x));
addOptional(ip, 'NA', 1.4, @(x)isnumeric(x));
addOptional(ip, 'ni', 1.5, @(x)isnumeric(x));
addOptional(ip, 'ns', 1.41, @(x)isnumeric(x));
addOptional(ip, 'lambda', 610, @(x)isnumeric(x));
addOptional(ip, 'microscope', 'WideField', @(x)ischar(x));

parse(ip, varargin{:});
ip = ip.Results;

ip.ng = 1.525;
ip.tg = 120000;
ip.tgd = 120000;

% Scaling
scale = 'linear';

% Type of microscope used
if strcmp(ip.microscope,'WideField') || strcmp(ip.microscope,'SpinDisk')
    n = 0;
elseif strcmp(ip.microscope,'Confocal')
    n = 1;
else 
    error('no valid microscope selected');
end

% Initialize the subvolume
value = zeros(ip.n);
nx = ip.n(1);
ny = ip.n(2);
nz = ip.n(3);

% Define center of subvolume
x0 = (nx+1)/2;
y0 = (ny+1)/2;
z0 = (nz+1)/2;

% Optical axis should go through the center
xp = x0;
yp = y0;


% lambda = 450;
% zp = 2000;

% Microscope tube length
zd = 1000;
zds = zd;


% Further parameter that are not really needed
M = 40; % Magnification
A0 = 1; % Constant
ti0 = 150*1e3; % Design working distance

% Resolution 
resLateral = ip.res(1);
resAxial = ip.res(3);
k0 = 2*pi/ip.lambda;

% xc = ((1:nx)-x0-1)*resLateral;
% yc = ((1:ny)-y0-1)*resLateral;
% zc = ((1:nz)-z0-1)*resAxial;

% Additional factor (unused)
a = ip.NA*zds/M;

% h  = (k*a.^2*A0/(zds.^2))^2 * abs(quad(@opd,0,1))^2;

if strcmp(ip.microscope,'SpinDisk')
    hexPattern = Simulation.getHexaPattern(ip.n,resLateral);
end

maxRadius = round(sqrt((nx-x0).^2+(ny-y0).^2))+1;
rad = zeros(1,maxRadius);
oversampling = 2;

% matlabpool open

for z = 1:nz
%     ti = ti0 + (z-z0)*resAxial; % not true -> no -1?
    
    for k = 1:maxRadius*oversampling
        r(k) = (k-1)*resLateral/oversampling;
        rad(k) = (k-1)/oversampling;
        h(k) = abs(integral(@opd,0,1)).^(2*n+2);
    end
    
    for x = 1:nx
        for y = 1:ny
            xd = x;
            yd = y;
            
            rPixel = sqrt((xd-xp).^2+(yd-yp).^2);
            rIndex = floor(rPixel*oversampling)+1;
%             value(x,y,z) = interp1(r,h,rPixel*resLateral);
            value(x,y,z) = h(rIndex) + (h(rIndex+1)-h(rIndex))*(rPixel-rad(rIndex))*oversampling;
%             h(x,y,z) = (k*a.^2*A0/(zds.^2))^2 *t;          
        end
    end
    
    % Spinning disk setup
    if strcmp(ip.microscope,'SpinDisk')
        value(:,:,z) = value(:,:,z).*conv2(value(:,:,z),hexPattern,'same');
    end
    
  %Show progress
  fprintf([repmat(sprintf('\b'),1,1+numel(num2str(round(100*z/nz)))),...
    sprintf('%d%%%%',round(100*z/nz))]);

end
% matlabpool close

    function I = opd(rho)
%         w = ip.ns*ip.zp*sqrt(1-(ip.NA*rho/ip.ns).^2) + ip.ni*(ti-ti0)*sqrt(1-(ip.NA*rho/ip.ni).^2) + a.^2*(zds-zd)/(2*zds*zd)*rho.^2 ; % Aguet 2013
%           achtung, aguet 2013 geht nicht weil noch nicht an neues z0
%           angepasst!
%         w = ((z0-z)*resAxial*sqrt(ip.ni^2-(ip.NA*rho).^2)+ip.zp(sqrt(ip.ns^2-(ip.NA*rho).^2)-(ip.ni/ip.ns).^2*sqrt(ip.ni^2-(ip.NA*rho).^2))); % Aguet 2005
%         w = ip.ns*ip.zp*sqrt(1-(ip.NA*rho/ip.ns).^2)+k0*ip.ni*(z0-z+1)*resAxial*sqrt(1-(ip.NA*rho/ip.ni).^2); % Shit
        w = ip.ni*(z0-z)*resAxial*sqrt(1-(ip.NA*rho/ip.ni).^2)+ip.ns*ip.zp*(sqrt(1-(ip.NA*rho/ip.ns).^2)-(ip.ni/ip.ns).^2*sqrt(1-(ip.NA*rho/ip.ni).^2)); % Gibson Lanny 1991 (same results as aguet 2005)
%         w = ip.ns*ip.zp*sqrt(1-(ip.NA*rho/ip.ns).^2)-ip.zp*ip.ni*sqrt(1-(ip.NA*rho/ip.ni).^2);
%         w = ip.ns*ip.zp*(sqrt(1-(ip.NA*rho/ip.ns).^2)-(ip.ni/ip.ns).^2*sqrt(1-(ip.NA*rho/ip.ni).^2)); %
%         w = ip.ni*(z0-z)*resAxial*sqrt(1-(ip.NA*rho/ip.ni).^2)+ip.ns*ip.zp*(sqrt(1-(ip.NA*rho/ip.ns).^2)-(ip.ni/ip.ns).^2*sqrt(1-(ip.NA*rho/ip.ni).^2))+ ...
%             +ip.ng*(sqrt(1-(ip.NA*rho/ip.ng).^2)-(ip.ni/ip.ng).^2*sqrt(1-(ip.NA*rho/ip.ni).^2))*(ip.tg-ip.tgd); % Gibson Lanny 1991 (same results as aguet 2005) with glass
        w = k0*w;

        I = besselj(0,k0*ip.NA*r(k)*rho).*exp(1i*w).*rho;
    end



[maxVal, I] = max(value(:));

% linear scaling
out =value/maxVal;
if strcmp(scale,'log')
% Logarithmic scaling
    out = log(out);
end

fprintf('\n')

end