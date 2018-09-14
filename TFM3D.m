function  TFM3D()
%TFM3D is a framework for 3D TFM simulation or analysis of experimental data
%CL

clear;
warning off;

%*******************************************************************************
% Set parameters
%*******************************************************************************

% Mode
%-------------------------------------------------------------------------------
mode = 'i2u';
% 'simt' Full 3D TFM simulation                  FEA data --> images --> u --> T
% 'simu' Generate images, calc. displacements    FEA data --> images --> u
% 'simi' Generate images                         FEA data --> images
% 'expu' Get 3D TFM images, calc. displacements  experimental images --> u
% 'expt' Full experimental 3D TFM analyis        experimental images --> u --> T
% 'i2t'  Get images, calc. tractions                          images --> u --> T
% 'i2u'  Get images, calc. displacements                      images --> u
% 'u2t'  Get displacements, calc. tractions                              u --> T
                         
% Image generation parameters for 3D TFM simulation (for mode = 'sim*')
%---------------------------------------------------------  ----------------------
sim.folder = '5_Spheroid21_x2';         % Project folder in 'TFM 3D\FEAData'
sim.cbeads = [ 0.1 ];                   % Bead concentration [1/µm^3]       [1xn]
sim.rbead = [0.1];                     % Bead radius [µm]              [1xn]
sim.vxsize = [0.2,0.2,0.2];            % Voxel size [µm]                   [nx3]
sim.PSFparam = [1.2,1.33,1.33];        % PSF: [NA,ni,ns], [] = PSF off     [nx3]
sim.ndatasets = 1;                      % Number of datasets to be generated

% Parameters for analsysis based on experimental data (for mode = 'expt')
%-------------------------------------------------------------------------------
exp.folder = '2_Cell5';                % Project folder in 'TFM 3D\ExpData'
exp.trsh = 1;                          % Intensity threshold [0,1]*Imax
exp.rbead = 0.1;                       % Bead radius [µm]
exp.vxsize = [0.2841,0.2841,0.8];      % Voxel size [µm]                   [1x3]
exp.PSFparam = [1.1,1.33,1.33];        % PSF: [NA,ni,ns] if available      [1x3]

% Displacement calculation parameters
%-------------------------------------------------------------------------------
ualg = {'FIDVC'};      % Displacement calculation algorithms:      {1xn}
                              %   'PDVC'   particle-based DVC
                              %   'GDVC'   grid-based DVC
                              %   'known'  include simulated displacements
PDVC.umax = [4];              % PDVC: expected maximum displacement [µm]   [1xn]
PDVC.cbeads_win = [0.1];      % PDVC: bead concentration per search window [1xn]
GDVC.maxwinsize = [64,64,64]; % GDVC: initial search window size [px]      [nx3]
GDVC.minwinsize = [32,32,32]; % GDVC: minimum search window size [px]      [nx3]
GDVC.spacing = [24,24,24];       % GDVC: spacing of search windows [px]       [nx3]
FIDVC.sSize = [128 128 64]; %FIDVC Interrogation window size [px]      [nx3]
                                
itp = {'grid'}; % Displacement interpolation locations:      {1xn}
                              %   'bead'   detected neutral bead positions
                        	  %   'grid'   regular grid
                              %   'cell'   cell vertices
itp_spacing = [1];          % Grid spacing for itp = 'grid' [µm]         [1xn]
ROIsize = [NaN];          % Region of interest for itp [% cell size]   [1xn]
                              %   (NaN = full-field analysis)

% Traction calculation parameters
%-------------------------------------------------------------------------------
E = [5000];                   % Young's modulus of the substrate [Pa]
nu = [0.45];                  % Poisson's ratio of the substrate
Talg = {'Direct'};
                              % Traction calculation algorithms:           {1xn}
                              %   'Direct' Direct differentiation of u
                              %   'LsqFit' Differentiation of LSQ fitted u
                              %   'FEA'    Finite element analysis
                              %   'Green'  Solve u=GT with Green's function
Direct.spacing = [1];         % Direct: grid spacing [µm]                  [1xn]
LsqFit.spacing = [1];         % LsqFit: grid spacing [µm]                  [1xn]
LsqFit.winsize = [3];         % LsqFit: window size for grad(u) fit [-]    [1xn]
FEA.elemsize = [0.5,10];        % FEA: surface element size [cell,gel] [µm]  [nx2]  
Green.elemsize = [5,50];      % Green: surface element size [cell,gel] [µm][nx2]

%*******************************************************************************
% Do not edit!
%*******************************************************************************

% Determine the number of simulation series and datasets
if strcmp('sim',mode(1:3)) && length(mode) == 4 % For 3D TFM simulation
    % 1 series of datasets per specified bead density
    ndatasets = sim.ndatasets;
    nseries = size(sim.cbeads,2);
elseif strcmp('exp',mode(1:3)) % For 3D TFM analysis
    ndatasets = 1;
    nseries = 1;
elseif strcmp('2',mode(2)) % For calculations based on imported datasets
    % Select folder to load all existing datasets from
    folder_res = uigetdir(strcat(pwd,'\Results'),'Select Input Folder');
    fileslist = dir(strcat(folder_res,'\*dataset*.mat'));
    ndatasets = size(fileslist,1);
    nseries = 1;
else
    fprintf('Error: Invalid TFM3D mode.\n')
    return
end

% Group displacement field calculation parameters
uparam.ualg = ualg;
uparam.itp = itp;
if any(strcmp('grid',itp)) || any(strcmp('GDVC',ualg)) || any(strcmp('FIDVC',ualg))
    uparam.itp_spacing = itp_spacing;
end;
uparam.ROIsize = ROIsize;
if any(strcmp('PDVC',ualg)), uparam.PDVC = PDVC; end;
if any(strcmp('GDVC',ualg)), uparam.GDVC = GDVC; end;
if any(strcmp('FIDVC',ualg)), uparam.FIDVC = FIDVC; end;
if any(strcmp('known',ualg)), uparam.known.incl = true; end;

% Group traction field calculation parameters
Tparam.E = E;
Tparam.nu = nu;
Tparam.Talg = Talg;
if any(strcmp('Direct',Talg)), Tparam.Direct = Direct; end;
if any(strcmp('LsqFit',Talg)), Tparam.LsqFit = LsqFit; end;
if any(strcmp('FEA',Talg)), Tparam.FEA = FEA; end;
if any(strcmp('Green',Talg)), Tparam.Green = Green; end;

% Global results structure
serieseval = [];
simtime = 0;
disptime.PDVC = 0;
disptime.GDVC = 0;
disptime.FIDVC = 0;
tractime = 0;
evaltime = 0;

% Generate a new simulation series for each given bead concentration
for i = 1:nseries
    
    % Generate multiple datasets for statistical evaluation
    for j = 1:ndatasets
        
        tic
        
        % Generate or load 3D images
        if strcmp('sim',mode(1:3)) % For 3D TFM simulation
            
            % Generate image pairs with identical bead distribution bundled as a dataset
            fprintf('Series %d/%d, dataset %d/%d:\n',i,nseries,j,ndatasets)
            %----------------------------------------------------------------------------------------
            dataset = CreateBeadImageSet(sim.folder,sim.cbeads(i),sim.rbead,sim.vxsize,sim.PSFparam);
            %----------------------------------------------------------------------------------------
            datasetlbl{j} = sprintf('dataset%03d',j);
            imgpairlbl = fieldnames(dataset);
            nimgpairs = size(imgpairlbl,1);
            
        elseif strcmp('exp',mode(1:3)) % For 3D TFM analysis
            
            % Create dataset from experimental TFM image data
            fprintf('Load 3D image pair and cell geometry... ')
            datasetlbl = {'dataset001'};
            imgpairlbl = {'imgpair001'};
            nimgpairs = 1;
            name = exp.folder;
            %------------------------------------------------------------------------
            img0 = ReadTIFFStack(strcat(pwd,'\ExpData\',name,'\',name,'_img0.tif'));
            img1 = ReadTIFFStack(strcat(pwd,'\ExpData\',name,'\',name,'_img1.tif'));
            exp.fvcell = load(strcat(pwd,'\ExpData\',name,'\',name,'.mat'));
            %------------------------------------------------------------------------
            % Simple filtering
            img0(img0<exp.trsh*max(img0(:))) = 0;
            img1(img1<exp.trsh*max(img1(:))) = 0;
            dataset.(imgpairlbl{1}).img0 = img0;
            dataset.(imgpairlbl{1}).img1 = img1;
            % Store parameters
            range = size(img0).*exp.vxsize;
            exp.range = [-range'/2,range'/2]; % [0,0,0] at image center
            dataset.(imgpairlbl{1}).imgparam = struct('projectname',exp.folder,'rbead',exp.rbead,...
                'vxsize',exp.vxsize,'PSFparam',exp.PSFparam,'range',exp.range,'fvcell',exp.fvcell);
            fprintf('done\n')

        else % For calculations based on imported datasets
            
            % Load existing dataset from the selected folder
            filename = strcat(folder_res,'\',fileslist(j).name);
            fprintf('Loading dataset %d/%d: %s... ',j,ndatasets,filename)
            %-----------------------------
            dataset_temp = load(filename);
            %-----------------------------
            datasetlbl{j} = sprintf('dataset%03d',j);
            imgpairlbl = fieldnames(dataset_temp);
            nimgpairs = size(imgpairlbl,1);
            
            % Keep only specified input data
            for n = 1:nimgpairs
                dataset.(imgpairlbl{n}).imgparam = dataset_temp.(imgpairlbl{n}).imgparam;
                dataset.(imgpairlbl{n}).img0 = dataset_temp.(imgpairlbl{n}).img0;
                dataset.(imgpairlbl{n}).img1 = dataset_temp.(imgpairlbl{n}).img1;
                if strcmp('u2t',mode)
                    dataset.(imgpairlbl{n}).uparam = dataset_temp.(imgpairlbl{n}).uparam;
                    dataset.(imgpairlbl{n}).u = dataset_temp.(imgpairlbl{n}).u;
                end;
            end
            fprintf('done\n')
            
        end
            simtime = toc;
        % Calculation for each image pair of the dataset
        for k = 1:nimgpairs
%             tic
            if ~strcmp('simi',mode), fprintf('*** Image pair %d/%d ***\n',k,nimgpairs); end;
            
            % Calculate displacements and attach the results to the data structure
            if ~strcmp('simi',mode) && ~strcmp('u2t',mode)
                %------------------------------------------------
                [u, disptime] = Displacement(dataset.(imgpairlbl{k}),uparam);
                %------------------------------------------------
                dataset.(imgpairlbl{k}).uparam = uparam;
                dataset.(imgpairlbl{k}).u = u;
            end
%             disptime = toc;
            tic
            % Calculate traction stresses and attach the results to the data structure
            if strcmp('t',mode(end))
                %--------------------------------------------
                T = Traction(dataset.(imgpairlbl{k}),Tparam);
                %--------------------------------------------
                dataset.(imgpairlbl{k}).Tparam = Tparam;
                dataset.(imgpairlbl{k}).T = T;
            end
            tractime = toc;
        end
        tic
        % Add results of the current dataset to the respective results tables
        if strcmp('u',mode(end))
            serieseval = SeriesEval(serieseval,dataset,'u');
        elseif strcmp('t',mode(end))
            serieseval = SeriesEval(serieseval,dataset,'T');
        end
        
        % Save the dataset
        if j == 1 % Create a new folder for a new series of datasets
            foldername = strcat('Series-',datestr(now,'yymmddHHMMSS'));
            mkdir('.\Results',foldername);
        end
        filename = sprintf('.\\Results\\%s\\%s.mat',foldername,datasetlbl{j});
        save(filename,'-struct','dataset','-v7.3');
        fprintf('*** Dataset saved to %s%s ***\n',pwd,filename(2:end))
        evaltime = toc;
        % Write timeresults to file
        fileID = fopen(sprintf('.\\Results\\%s\\time_results.txt',foldername),'a');
%         fprintf(fileID,'%s %s %s %s %s\n','datasetlbl{j}','simtime','disptime','tractime','evaltime');
        fprintf(fileID,'%s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\r\n',datasetlbl{j},simtime,disptime.PDVC,disptime.GDVC,disptime.FIDVC,tractime,evaltime);
        fclose(fileID);
    end
    
    % Save series evaluation data
    if ~isempty(serieseval)
        filename_eval = sprintf('.\\Results\\%s\\serieseval.mat',foldername);
        save(filename_eval,'-struct','serieseval');
        fprintf('*** Series evaluation saved to %s%s ***\n\n',pwd,filename_eval(2:end))
    end
    
end

if strcmp('u',mode(end))
    ResultsViewer(foldername,'u')
elseif strcmp('t',mode(end))
    ResultsViewer(foldername,'T')
end

warning on;

end
