function CompareFEAResults()
%COMPAREFEARESULTS determines a modified mean absolute percentage error (MAPE)
%                  between two ANSYS output data sets containing displacements
%                  and normal stresses in x,y,z
%CL

% Set parameters
%-----------------------------------------------------------------------------
mode1 = 'auto'; % 'manual' or 'auto' mode of FEA results export 1
mode2 = 'auto'; % 'manual' or 'auto' mode of FEA results export 2
spacing = 1;    % Grid spacing for interpolation
ROI = [-2,2;-2,2;12,16]; % Region of interest [xmin,xmax;ymin,ymax;zmin,zmax]
%-----------------------------------------------------------------------------

warning off

fprintf('Spacing = %f, ROI = [%.1f,%.1f;%.1f,%.1f;%.1f,%.1f]\n',spacing,ROI');

% Get FEA data folder 1
folder1 = uigetdir(strcat(pwd,'\FEAData'));
[~,projectname1,~] = fileparts(folder1);
% Get FEA data folder 2
folder2 = uigetdir(strcat(pwd,'\FEAData'));
[~,projectname2,~] = fileparts(folder2);

% Load data
fprintf('Loading results from %s\\, mode: %s ',folder1,mode1)
[fuFEA1,fsFEA1,range1] = ImportFEAResults(strcat(folder1,'\',projectname1,'.wbpj'),mode1);
fprintf('\nLoading results from %s\\, mode: %s ',folder2,mode2)
[fuFEA2,fsFEA2,range2] = ImportFEAResults(strcat(folder2,'\',projectname2,'.wbpj'),mode2);
fprintf('\nRange 1 = [%.1f,%.1f;%.1f,%.1f;%.1f,%.1f]\n',range1')
fprintf('Range 2 = [%.1f,%.1f;%.1f,%.1f;%.1f,%.1f]\n',range2')

% Interpolation on cubic grid with a range fully inside the range of FE data
[x,y,z] = meshgrid(ROI(1,1):spacing:ROI(1,2),ROI(2,1):spacing:ROI(2,2),ROI(3,1):spacing:ROI(3,2));
nodes = [x(:),y(:),z(:)];
dim = {'x','y','z'};

for i = 1:3 % Interpolate displacement results
    u1(:,i) = fuFEA1.(dim{i})(nodes);
    u2(:,i) = fuFEA2.(dim{i})(nodes);
end
uMAPE = mean(100*norm(u1-u2)/max(norm(u1),norm(u2)));
fprintf('uMAPE = %f\n',uMAPE)

for j = 1:3 % Interpolate stress results
    s1(:,j) = fsFEA1.(dim{j})(nodes);
    s2(:,j) = fsFEA2.(dim{j})(nodes);
end
sMAPE = mean(100*norm(s1-s2)/max(norm(s1),norm(s2)));
fprintf('sMAPE = %f\n',sMAPE)

warning on

end