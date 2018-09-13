function [ spos, scalc ] = FEAStressCalc( upos, ucalc, E, nu, elemsize_cell, elemsize_sub, range, fvcell )
%FEASTRESSCALC      calculates the stress field based on known displacements
%                   using a finite element analysis in ANSYS Mechanical APDL.
%Input:
%  <upos>           positions where displacements are known [µm]
%  <ucalc>          known displacements [µm]
%  <E>              Young's modulus of the substrate [Pa]
%  <nu>             Poisson's ratio of the substrate
%  <elemsize_cell>  cell surface element size [µm]
%  <elemsize_sub>   substrate surface element size [µm]
%  <range>          nodal range [µm]
%  <fvcell>         cell model patch object [µm]
%Output:
%  <spos>           positions where stress tensors are calculated [µm]
%  <s>              stress tensors in Nye notation [sx,sy,sz,sxy,sxz,syz] [Pa]
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

% Write displacements to input file for ANSYS
npoint = 0:size(upos,1);
utable = table(npoint',[1,2,3,4,5,6;upos,ucalc]);
writetable(utable,strcat(work_dir,'\u.txt'),'WriteVariableNames',false);

% Write coordinate range to input file for ANSYS
trange = table([0;1;2;3],[1,2;range]);
writetable(trange,strcat(work_dir,'\range.txt'),'WriteVariableNames',false);

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
fid = fopen(strcat(pwd,'\Traction\FEA\FEA2StressCalcTemplate.cmd'));
n = 0;
while ~feof(fid)
    n=n+1;
    tline{n,1} = fgetl(fid);
end
fclose(fid);

% Create ANSYS input file
inputfile = strcat(work_dir,'\FEAStressCalc.cmd');
fid = fopen(inputfile,'wt');
for i = 1:size(tline,1)
    if strfind(tline{i},'/CWD')
        % Write working directory
        fprintf(fid,'/CWD,''%s''\n',work_dir);
    elseif strfind(tline{i},'~SATIN')
        % Write path to cell model
        fprintf(fid,'~SATIN,''%s'',''%s'',''%s'',SOLIDS,0\n',name_cell,ext_cell,path_cell);
%     elseif strfind(tline{i},'BLOCK')
%         % Write command to create the substrate block
%         fprintf(fid,'BLOCK,%f,%f,%f,%f,%f,%f\n',range(1,:),range(2,:),range(3,:));
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

% Wait for analysis to be completed, i.e. until results file is written
outputfile = strcat(work_dir,'\stress.txt');
fid = -1;
while fid == -1
    fid = fopen(outputfile);
end
fclose(fid);

% Delete batch file
delete(batch_file);

% Read ANSYS output file
pause(5);
data = csvread(outputfile,1,1);
spos = data(:,2:4); % Nodal coordinates of the finite element model
scalc = data(:,5:10); % Stress tensors at nodes of the finite element model


end

