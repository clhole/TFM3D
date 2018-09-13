function [ dataset ] = CreateBeadImageSet( projectname, cbeads, rbead, vxsize, PSFparam )
%CREATEBEADIMAGESET creates 3D bead image pairs based on the displacement
%                   results of an ANSYS FEM simulation of an embedded cell.
%                   Multiple inputs for bead radius and voxel size are
%                   allowed and result in multiple image pairs of the
%                   same bead distribution that are bundled as a dataset.
%Input:
%  <projectname>    name of the folder in 'TFM 3D\FEAData\' containing the
%                   same-named .stl cell geometry and ANSYS project files
%  <cbeads>         bead concentration [beads/µm^3]
%  <rbead>          bead radius [µm]                                       [1xn]
%  <vxsize>         voxel size [dx,dy,dz] [µm]                             [nx3]
%  <PSFparam>       Point spread function parameters: [NA,ni,ns]           [nx3]
%                      NA         numerical aperture
%                      ni         refractive index of the immersion liquid
%                      ns         refractive index of the substrate
%                   (Set PSFparam = [] to turn off the PSF.)
%Output:
%  <dataset>        structure containing the data of each image pair:
%                   .imgparam
%                      .cbeads    bead concentration [beads/µm^3]
%                      .rbead     bead radius [µm]
%                      .vxsize    voxel size [µm]
%                      .range     nodal range [µm]
%                      .fvcell    cell model patch object [µm]
%                      .fuFEA     FEA displacement interpolant [µm]
%                      .fTFEA     FEA stress interpolant [Pa]
%                      .beadpos0  given bead positions [µm] at t=t0 
%                      .beadpos1  given bead positions [µm] at t=t1
%                   .img0         bead image at t=t0 (with cell)
%                   .img1         bead image at t=t1 (without cell)
%CL

fprintf('Importing FEA data')

% Import cell surface patch from .mat file
file_cell = strcat(pwd,'\FEAData\',projectname,'\',projectname,'.mat');
fvcell = load(file_cell);

% Import FEA results from ANSYS Workbench (as interpolants)
file_wb = strcat(pwd,'\FEAData\',projectname,'\',projectname,'.wbpj');
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
fprintf(' %d bead positions generated.\n',nbeads)

nimgpairs = 0;

% Repeat for multiple voxel sizes
for i = 1:size(vxsize,1)
        
    % Repeat for multiple bead radii
    for j = 1:size(rbead,2)
        
        % Repeat for each set of PSF parameters
        for k = 1:size(PSFparam,1)
        
            nimgpairs = nimgpairs + 1;
            name = sprintf('imgpair%03d',nimgpairs);
            
            % Create image pair
            PSFparamset = PSFparam(k,:);
            img0 = Create3DImage(vxsize(i,:),range,rbead(j),beadpos0,PSFparamset);
            img1 = Create3DImage(vxsize(i,:),range,rbead(j),beadpos1,PSFparamset);
            
            % Create dataset
            dataset.(name).imgparam = struct('projectname',projectname,...
                'cbeads',cbeads,'rbead',rbead(j),'vxsize',vxsize(i,:),...
                'PSFparam',PSFparamset,'range',range,'fvcell',fvcell,'fuFEA',fuFEA,...
                'fTFEA',fTFEA,'frFEA',frFEA,'beadpos0',beadpos0,'beadpos1',beadpos1);
            dataset.(name).img0 = img0;
            dataset.(name).img1 = img1;
            
        end
        
        if isempty(PSFparam) % PSF off
            
            nimgpairs = nimgpairs + 1;
            name = sprintf('imgpair%03d',nimgpairs);
            
            % Create image pair
            img0 = Create3DImage(vxsize(i,:),range,rbead(j),beadpos0);
            img1 = Create3DImage(vxsize(i,:),range,rbead(j),beadpos1);
            
            % Create dataset
            dataset.(name).imgparam = struct('projectname',projectname,...
                'cbeads',cbeads,'rbead',rbead(j),'vxsize',vxsize(i,:),...
                'PSFparam',[],'range',range,'fvcell',fvcell,'fuFEA',fuFEA,...
                'fTFEA',fTFEA,'frFEA',frFEA,'beadpos0',beadpos0,'beadpos1',beadpos1);
            dataset.(name).img0 = img0;
            dataset.(name).img1 = img1;
            
        end
       
    end
    
end

fprintf('Dataset containing %d image pair(s) created.\n',nimgpairs)

end