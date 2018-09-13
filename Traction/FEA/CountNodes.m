function [ nnodes ] = CountNodes( filename, range, elemsize_min, elemsize_max )
%COUNTNODES     	returns the total number of corner nodes of the FEA
%                   traction calculation setup in ANSYS Mechanical APDL.
%Input:
%   <filename>  	geometry file in .sat format
%	<range>     	coordinate range [xmin,xmax;ymin,ymax;zmin,zmax] of a
%               	cuboid from which the .sat geometry gets subtracted    [1x3] 
%	<elemsize_min>  minimum element size [µm] (substrate edges)
%	<elemsize_max>  maximum element size [µm] (inner substrate boundary)
%Output:
%	<nnodes>        number of (corner) nodes
%CL

warning off

% Set parameters
ext = 'sat';                                 % File extension
ANSYS_lic = 'aa_r';                          % ANSYS: Use academic research license
work_dir = strcat(pwd,'\Traction\FEA\work'); % Working directory

% Get filename
[path_sat,name_sat,~] = fileparts(filename);

% Delete working directory and create a new one
fclose all;
try
    rmdir(work_dir,'s');
catch
end
mkdir(work_dir);

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
fid = fopen(strcat(pwd,'\Traction\FEA\CountNodesTemplate.cmd'));
    n = 0;
    while ~feof(fid)
        n=n+1;
        tline{n,1} = fgetl(fid);
    end
fclose(fid);

% Create ANSYS input file
inputfile = strcat(work_dir,'\CountNodes.cmd');
fid = fopen(inputfile,'wt');
for i = 1:size(tline,1)
    if strfind(tline{i},'/CWD')
        % Write working directory
        fprintf(fid,'/CWD,''%s''\n',work_dir);
    elseif strfind(tline{i},'~SATIN')
        % Write path to cell model
        fprintf(fid,'~SATIN,''%s'',''%s'',''%s'',SOLIDS,0\n',name_sat,ext,path_sat);
    elseif strfind(tline{i},'BLOCK')
        % Write command to create the substrate block
        fprintf(fid,'BLOCK,%f,%f,%f,%f,%f,%f\n',range(1,:),range(2,:),range(3,:));
    elseif strfind(tline{i},'VGEN')
        % Write command to shift the cell to the substrate center
        fprintf(fid,'VGEN,,1,,,%f,%f,%f,,,1\n',range(:,1)+(range(:,2)-range(:,1))/2);
    elseif strfind(tline{i},'LESIZE,edges')
        % Write mesh spacing for substrate edges
        fprintf(fid,'LESIZE,edges,%.1f\n',elemsize_max);
    elseif strfind(tline{i},'ESIZE')
        % Write mesh spacing for inner substrate boundary
        fprintf(fid,'ESIZE,%.1f\n',elemsize_min);
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
fprintf(fid,'"%s" -p %s -dir "%s" -j "CountNodes" -b -i "%s" -o "%s\\CountNodes.out"\n',...
    ANSYS_path,ANSYS_lic,work_dir,inputfile,work_dir);
fprintf(fid,'EXIT');
fclose(fid);

% Run batch file via Windows start menu as a workaround to prevent ANSYS from
% producing a solver stack overflow when being directly called from within MATLAB.
RunBat(batch_file);

% Determine if the test was successfull: Check for existence of output file
outputfile = strcat(work_dir,'\CountNodes.txt');
fid = -1;
tic
while fid == -1
    fid = fopen(outputfile);
    if fid ~= -1 || toc > 15 % File found or reached time limit
        break
    end
end
if fid == -1
    fprintf('CountNodes: Test failed. Read log files in TFM 3D\\Traction\\FEA\\work.\n')
else
    data = textscan(fid,'Nodes,%d');
    nnodes = data{:};
    fclose(fid);
end

% Delete batch file
delete(batch_file);

warning on

end

