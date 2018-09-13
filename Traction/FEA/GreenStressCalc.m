function [ Tpos, Tcalc ] = GreenStressCalc( upos, ucalc, E, nu, elemsize_cell, elemsize_sub, fvcell )
%GREENSTRESSCALC    calculates the stress field by generating a discretized
%                   Green's function G which relates tractions T on the
%                   cell surface to displacements u within the surrounding
%                   substrate using ANSYS Mechanical APDL.
%                   Method to generate G: Wesley Legant, 2010
%                   Optimization to solve u=GT: PC Hansen, 2003                
%Input:
%  <upos>           positions where displacements are known [µm]
%  <ucalc>          known displacements [µm]
%  <E>              Young's modulus of the substrate [Pa]
%  <nu>             Poisson's ratio of the substrate
%  <elemsize_cell>  cell surface element size [µm]
%  <filename_cell>  cell surface patch file (.mat)
%  <fvcell>         cell model patch object [µm]
%Output:
%  <Tpos>           positions where the traction stresses are calculated [µm]
%  <Tcalc>          calculated traction stresses [Pa]
%CL

% Set parameters
work_dir = strcat(pwd,'\Traction\FEA\work'); % Working directory
ANSYS_lic = 'aa_r'; % ANSYS: Use academic research license
  
% Delete working directory and create a new one
fclose all;
try
    rmdir(work_dir,'s');
catch
end
mkdir(work_dir);

% Generate .sat file for input to ANSYS Mechanical APDL
path_cell = work_dir;
name_cell = 'cellmodel';
ext_cell = 'sat';
filename_sat = strcat(path_cell,'\',name_cell,'.',ext_cell);
fv2sat(fvcell,filename_sat)
nfaces = size(fvcell.faces,1);

% Write query point locations to input file for ANSYS
npoint = 0:size(upos,1);
upostable = table(npoint',[1,2,3;upos]);
writetable(upostable,strcat(work_dir,'\upos.txt'),'WriteVariableNames',false);

% Get path of the newest ANSYS version
for i = 1:0.5:30
   directory = getenv(sprintf('ANSYS%d_DIR',i*10));
   if ~strcmp(directory,'')
       ANSYS_dir = directory;
       ANSYS_vers = i;
   end
end
sys = getenv('ANSYS_SYSDIR');
ANSYS_path = sprintf('%s\\bin\\%s\\ANSYS%d.exe',ANSYS_dir,sys,ANSYS_vers*10);

% Read ANSYS input template file
fid = fopen(strcat(pwd,'\Traction\FEA\GreenStressCalcTemplate.cmd'));
n = 0;
while ~feof(fid)
    n=n+1;
    tline{n,1} = fgetl(fid);
end
fclose(fid);

% Create ANSYS input file
inputfile = strcat(work_dir,'\GreenStressCalc.cmd');
fid = fopen(inputfile,'wt');
for i = 1:size(tline,1)
    if strfind(tline{i},'/CWD')
        % Write working directory
        fprintf(fid,'/CWD,''%s''\n',work_dir);
    elseif strfind(tline{i},'~SATIN')
        % Write path to cell model
        fprintf(fid,'~SATIN,''%s'',''%s'',''%s'',SOLIDS,0\n',name_cell,ext_cell,path_cell);
    elseif strfind(tline{i},'MP,EX')
        % Write Young's modulus in [Pa]
        fprintf(fid,'MP,EX,1,%.3f\n',E);
    elseif strfind(tline{i},'MP,PRXY')
        % Write Poisson's ratio 
        fprintf(fid,'MP,PRXY,1,%.3f\n',nu);
    elseif strfind(tline{i},'LESIZE,edges')
        % Write mesh spacing for substrate edges
        fprintf(fid,'LESIZE,edges,%.1f\n',elemsize_sub);
    elseif strfind(tline{i},'ESIZE')
        % Write mesh spacing for inner substrate boundary
        fprintf(fid,'ESIZE,%.1f\n',elemsize_cell);
    elseif strfind(tline{i},'***nareas-1***')
        fprintf(fid,'(%d(E16.9'','')E16.9)\n',3*nfaces-1);
    else
        fprintf(fid,strcat(tline{i},'\n'));
    end
end
fclose(fid);

% Create batch file for running ANSYS with the specified parameters
startmenu_dir = strcat(getenv('APPDATA'),'\Microsoft\Windows\Start Menu\Programs');
batch_file = strcat(startmenu_dir,'\RUNANSYS.bat');
fid = fopen(batch_file,'wt');
fprintf(fid,'@ECHO off\n');
fprintf(fid,'"%s" -p %s -dir "%s" -j "TFM3D" -b -i "%s" -o "%s\\TFM3D.out"\n',...
    ANSYS_path,ANSYS_lic,work_dir,inputfile,work_dir);
fprintf(fid,'EXIT');
fclose(fid);

% Run batch file via Windows start menu as a workaround to prevent ANSYS from
% producing a solver stack overflow when being directly called from within MATLAB.
RunBatNew(batch_file);

% Wait for the analysis to be completed, then import the results
%**************************************************************************

progressfile_loads = strcat(work_dir,'\progress.tmp');
progressfile_solve = strcat(work_dir,'\TFM3D.mntr');
green_file = strcat(work_dir,'\green.txt');
fcent_file = strcat(work_dir,'\fcent.txt');

% Wait for load steps to be defined for each face
progress = 0;
clearmsg = '';
while progress<100
    fid = -1;
    while fid == -1 % Wait for progress file to be written
        fid = fopen(progressfile_loads);
    end
    data = textscan(fid,'Face,%f,%f');
    fclose(fid);
    area = data{1};
    nareas = data{2};
    progress = area/nareas*100; 
    % Show progress
    msg = sprintf('Preparing faces... %d/%d',area,nareas);
    fprintf([clearmsg,msg])
    clearmsg = repmat(sprintf('\b'),1,numel(msg));
end
fprintf('\n')

% Wait for solution to be completed
progress = 0;
clearmsg = '';
while progress<100
    fid = -1;
    n = 0;
    while fid == -1 % Wait for .mntr file to be written
        fid = fopen(progressfile_solve);
    end
    while ~feof(fid) % Read complete .mntr file
        n=n+1;
        txtline{n,1} = fgetl(fid);
    end
    fclose(fid);
    ls = sscanf(txtline{end,1},'%d%'); % Read load step number (last line)
    if isempty(ls), ls = 0; end;
    progress = ls/(3*nareas)*100;
    % Show progress
    msg = sprintf('Solving load step... %d/%d',ls,3*nareas);
    fprintf([clearmsg,msg])
    clearmsg = repmat(sprintf('\b'),1,numel(msg));
end
fprintf('\n')

% Wait for results file to be written
fid = -1;
while fid == -1
    fid = fopen(green_file);
end
fclose(fid);

% Delete batch file
delete(batch_file);

% Import discrete Green's function from FEA results
pause(5);
G = csvread(green_file);

% Import face centroid locations from FEA results
Tpos = csvread(fcent_file);
% Check if face numbering was kept during FEA
poserr = Tpos-fvcell.UserData.facecent;
if max(poserr)>0.01
    fprintf('WARNING: Check face numbering. Results might be incorrect!')
end

% Solve u=GT
%**************************************************************************

% Create displacement vector [ux1;uy1;uz1;ux2;uy2;uz2;...]
u = reshape(ucalc',[],1);

% Singular value decomposition of G
[U,s,V] = csvd(G);

% Tikhonov regularization with L-curve criterion
lambda_v = logspace(-3,1,30);
[lambda,~,~,~] = l_curve(U,s,u,'Tikh',lambda_v,V);

% Solve inverse problem 
[F,~,~] = tikhonov(U,s,V,u,lambda);

% Reshape traction force matrix
Tcalc = transpose(reshape(F,3,[]));

end

