function sat2ansTest( varargin )
%SAT2ANSTEST    checks if a .sat geometry file can be sucessfully imported
%               and used for boolean operations in ANSYS Mechanical APDL
%               (especially to check for memory-related SIG$SEGV errors).         
%Input:
%   <filename>  (optional) geometry file in .sat format
%CL

warning off

% Get file
if ~isempty(varargin)
    [path_sat,name_sat,~] = fileparts(varargin{1});
else
    [filename,path_sat] = uigetfile('*.sat');
    [~,name_sat,~] = fileparts(filename);
end

% Set parameters
ext = 'sat';                                 % File extension
ANSYS_lic = 'aa_r';                          % ANSYS: Use academic research license
work_dir = strcat(pwd,'\Traction\FEA\work'); % Working directory

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
fid = fopen(strcat(pwd,'\Traction\FEA\sat2ansTestTemplate.cmd'));
    n = 0;
    while ~feof(fid)
        n=n+1;
        tline{n,1} = fgetl(fid);
    end
fclose(fid);

% Create ANSYS input file
inputfile = strcat(work_dir,'\sat2ansTest.cmd');
fid = fopen(inputfile,'wt');
for i = 1:size(tline,1)
    if strfind(tline{i},'/CWD')
        % Write working directory
        fprintf(fid,'/CWD,''%s''\n',work_dir);
    elseif strfind(tline{i},'~SATIN')
        % Write path to cell model
        fprintf(fid,'~SATIN,''%s'',''%s'',''%s'',SOLIDS,0\n',name_sat,ext,path_sat);
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
fprintf(fid,'"%s" -p %s -dir "%s" -j "sat2ansTest" -b -i "%s" -o "%s\\sat2ansTest.out"\n',...
    ANSYS_path,ANSYS_lic,work_dir,inputfile,work_dir);
fprintf(fid,'EXIT');
fclose(fid);

% Run batch file via Windows start menu as a workaround to prevent ANSYS from
% producing a solver stack overflow when being directly called from within MATLAB.
RunBat(batch_file);

% Determine if the test was successfull: Check for existence of output file
outputfile = strcat(work_dir,'\sat2ansTest.txt');
fid = -1;
tic
while fid == -1
    fid = fopen(outputfile);
    if fid ~= -1 || toc > 15 % File found or reached time limit
        break
    end
end
if fid == -1
    fprintf('sat2ansTest: Test failed. Read log files in ''TFM 3D\\Traction\\FEA\\work''.\n\n')
else
    fclose(fid);
    fprintf('sat2ansTest: Test successfully completed.\n\n')
end

% Delete batch file
delete(batch_file);

warning on

end

