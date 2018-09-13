function [ fuFEA, fsFEA, frFEA, range ] = ImportFEAResults( filename, varargin )
%IMPORTFEARESULTS imports displacement and stress results from output
%                 text files that were created with ANSYS Workbench
%Input:
%   <filename>    filename of the ANSYS Workbench project (.wbpj)
%	<mode>        (optional) mode of ANSYS output format:
%                 'auto'    import macro generated results (default)
%                 'manual'  import manually exported results named as
%                           <filename_ux/uy/uz/sx/sy/sz/sxy/sxz/syz> (.txt)
%Output:
%	<fuFEA>       displacement field interpolants [x,y,z,ux,uy,uz] [µm]
%	<fsFEA>       stress field interpolants [x,y,z,Tx,Ty,Tz,Txy,Txz,Tyz] [Pa]
%	<range>       nodal coordinate range [x1,x2;y1,y2;z1,z2] [µm]
%CL

% Get path and name
[path_wb,name_wb,~] = fileparts(filename);

% Define mode
mode = 'auto';
if ~isempty(varargin)
    mode = varargin{1};
end

% Define dimensions
dim_u = {'x','y','z'};
dim_s = {'x','y','z','xy','xz','yz'};

if strcmp('auto',mode)
    
    % Import displacement results (ignore first row and column)
    file_u = strcat(path_wb,'\',name_wb,'_files\dp0\SYS\MECH\Results_u.txt');
    u = csvread(file_u,1,1);
    upos = u(:,2:4);
    uval = u(:,5:7);
    % Create interpolants
    for i = 1:size(dim_u,2)
        fuFEA.(dim_u{i}) = scatteredInterpolant(upos,uval(:,i));
        fprintf('.');
    end

    % Import stress results (ignore first row and column)
    file_s = strcat(path_wb,'\',name_wb,'_files\dp0\SYS\MECH\Results_s.txt');
    s = csvread(file_s,1,1);
    spos = s(:,2:4);
    sval = 10^6 * s(:,5:10); % Conversion from [MPa] to [Pa]
    % Create interpolants
    for i = 1:size(dim_s,2)
        fsFEA.(dim_s{i}) = scatteredInterpolant(spos,sval(:,i));
        fprintf('.');
    end
    
    % Import displacement results (ignore first row and column)
    % First File
    file_f1 = strcat(path_wb,'\',name_wb,'_files\dp0\SYS\MECH\input_f1.txt');
    res1 = csvread(file_f1,1,1);
    rpos1 = res1(:,2:4);
    rval1 = res1(:,5:7);
    for i = 1:size(dim_u,2)
        frFEA(1).(dim_u{i}) = scatteredInterpolant(rpos1,rval1(:,i));
        fprintf('.');
    end
    % Second File
    file_f2 = strcat(path_wb,'\',name_wb,'_files\dp0\SYS\MECH\input_f2.txt');
    res2 = csvread(file_f2,1,1);
    rpos2 = res2(:,2:4);
    rval2 = res2(:,5:7);
    for i = 1:size(dim_u,2)
        frFEA(2).(dim_u{i}) = scatteredInterpolant(rpos2,rval2(:,i));
        fprintf('.');
    end
    
    % Determine the range of nodal coordinates
    range(:,1) = min(u(:,2:4))';
    range(:,2) = max(u(:,2:4))';
    
else
    
    % Import displacement results
    for i = 1:3
        
        % Get the filename and open the corresponding file
        file_u = strcat(path_wb,'\',name_wb,'_u',dim_u{i},'.txt');
        fid = fopen(file_u);

        % Read header
        header = fgetl(fid);
        delimiter1 = strfind(header,'(');
        delimiter2 = strfind(header,')');

        % Determine the scaling factor for conversion to [µm]
        unit_length = header(delimiter1(1)+1:delimiter2(1)-1);
        switch unit_length
            case 'km', scf_length = 10^9;
            case 'm', scf_length = 10^6;
            case 'mm', scf_length = 10^3;
            case 'µm', scf_length = 1;
            otherwise, fprintf('FEA data error: Invalid unit of length.');
        end
        
        % Read data (consider that node numbers are optional output)
        data = textscan(fid,'%f %f %f %f %f');
        upos = scf_length*cell2mat(data(:,end-3:end-1));
        sval = scf_length*cell2mat(data(:,end));
        fuFEA.(dim_u{i}) = scatteredInterpolant(upos,sval);
        range_min(:,i) = min(upos);
        range_max(:,i) = max(upos);
        
        fclose(fid);
        fprintf('.');
        
    end
    
    % Import stress results
    for i = 1:6
        
        % Get the filename and open the corresponding file
        file_s = strcat(path_wb,'\',name_wb,'_s',dim_s{i},'.txt');
        fid = fopen(file_s);

        % Read header
        header = fgetl(fid);
        delimiter1 = strfind(header,'(');
        delimiter2 = strfind(header,')');

        % Determine the scaling factor for conversion to [µm]
        unit_length = header(delimiter1(1)+1:delimiter2(1)-1);
        switch unit_length
            case 'km', scf_length = 10^9;
            case 'm', scf_length = 10^6;
            case 'mm', scf_length = 10^3;
            case 'µm', scf_length = 1;
            otherwise, fprintf('FEA data error: Invalid unit of length.');
        end
        
        % Determine the scaling factor for conversion to [Pa]
        unit_stress = header(delimiter1(end)+1:delimiter2(end)-1);
        switch unit_stress
            case 'MPa', scf_stress = 10^6;
            case 'Pa', scf_stress = 1;
            otherwise, fprintf('FEA data error: Invalid unit of stress.');
        end
        
        % Read data (consider that node numbers are optional output)
        data = textscan(fid,'%f %f %f %f %f');
        spos = scf_length*cell2mat(data(:,end-3:end-1));
        sval = scf_stress*cell2mat(data(:,end));
        fsFEA.(dim_s{i}) = scatteredInterpolant(spos,sval);

        fclose(fid);
        fprintf('.');
        
    end
    
    % Determine the range of nodal coordinates
    range(:,1) = max(range_min,[],2)';
    range(:,2) = min(range_max,[],2)';

end

end

