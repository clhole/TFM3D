function [ spos, s ] = DirectStressCalc( upos, ucalc, E, nu, spacing )
%DIRECTSTRESSCALC directly differentiates the displacement field to obtain
%                 material strains which are then used to compute stresses
%Input:
%	<upos>        positions where the displacements are known [µm]
%	<ucalc>       known displacements [µm]
%	<E>           Young's modulus of the substrate [Pa]
%	<nu>          Poisson's ratio of the substrate
%   <spacing>     grid spacing for grad(u) calculation [µm]
%Output:
%   <spos>        positions where stress tensors are calculated [µm]
%	<s>           stress tensors in Nye notation [sx,sy,sz,sxy,sxz,syz] [Pa]
%CL

% Determine Lamé constants (homogeneous isotropic linear elastic material)
lambda = nu*E/((1+nu)*(1-2*nu));
% lambda = 0;
mu = E/(2*(1+nu));

% Prepare grid
range = [min(upos);max(upos)]';
xgv = ceil(range(1,1)*10)/10:spacing:floor(range(1,2)*10)/10;
ygv = ceil(range(2,1)*10)/10:spacing:floor(range(2,2)*10)/10;
zgv = ceil(range(3,1)*10)/10:spacing:floor(range(3,2)*10)/10;
[x,y,z] = meshgrid(xgv,ygv,zgv);
% Interpolate displacements
for d = 1:3
    fu = scatteredInterpolant(upos,ucalc(:,d));
    u{d} = fu(x,y,z);
end
nnodes = size(x(:),1);

% Smoothing
ux = smooth3(u{1},'box',3);
uy = smooth3(u{2},'box',3);
uz = smooth3(u{3},'box',3);
% ux = u{1};
% uy = u{2};
% uz = u{3};

% Differentiate displacement field
fprintf('\nDisplacement gradient calculation')
[duxdx,duxdy,duxdz] = gradient(ux,spacing);
fprintf('.')
[duydx,duydy,duydz] = gradient(uy,spacing);
fprintf('.')
[duzdx,duzdy,duzdz] = gradient(uz,spacing);
fprintf('.\n')

s = zeros(nnodes,6);
spos = zeros(nnodes,3);
n = 0;
clearmsg = '';

for i=1:size(x,1)
    for j=1:size(x,2)
        for k=1:size(x,3)
            
            n = n+1;
            
            % Calculate displacement gradient Du
            Du = [duxdx(i,j,k),duxdy(i,j,k),duxdz(i,j,k);
                  duydx(i,j,k),duydy(i,j,k),duydz(i,j,k);
                  duzdx(i,j,k),duzdy(i,j,k),duzdz(i,j,k)];
            
            % Calculate strain tensor
            epsilon = 1/2*(Du + Du');
            
            % Calculate stress tensor
            sigma = lambda*trace(epsilon)*eye(3) + 2*mu*epsilon;
            
            % Store stress tensor in Nye notation [sx,sy,sz,sxy,sxz,syz]
            s(n,:) = [sigma(1,1),sigma(2,2),sigma(3,3),sigma(1,2),sigma(1,3),sigma(2,3)];
         
            % Store grid node coordinates
            spos(n,:) = [x(i,j,k),y(i,j,k),z(i,j,k)];
            
        end
    end
    
    msg = sprintf('Stress field calculation... %.0f%%',i/size(x,1)*100);
    fprintf('%s',[clearmsg,msg])
    clearmsg = repmat(sprintf('\b'),1,length(msg));
    
end

fprintf('\n')

end

