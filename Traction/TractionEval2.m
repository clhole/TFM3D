function [ eval ] = TractionEval2( spos, scalc, imgparam, mode )
%TRACTIONEVAL   evaluates the stress field in the substrate (maximum
%               principal stress) and on the cell surface ('tractions')
%Input:
%   <spos>      positions where the stress field was calculated [µm]
%   <scalc>     stress tensors in Nye notation [sx,sy,sz,sxy,sxz,syz]
%               OR traction vectors [Tx,Ty,Tz]
%   <imgparam>  structure containing the following fields:
%               .fvcell>  cell model patch object incl. face centroids [µm]
%               .fTFEA>   (optional) given FEA stress interpolant [Pa]
%Output:
%   <eval>      table containing the evaluated traction stresses
%CL

% For simulated data, get given FEA stress interpolant for comparison
if isfield(imgparam,'fTFEA') 
    sim = true;
    fTFEA = imgparam.fTFEA;
else % For experimental data
    sim = false;
end

% Import FEA input
if isfield(imgparam,'frFEA')
    frFEA = imgparam.frFEA;
    inp_T = true;
end

switch mode
    case 'FDM'
        sgn = 1;
    case 'LSQ'
        sgn = 1;
    case 'FEA'
        sgn = -1;
    case 'GRE'
end

% Prepare stress field interpolants
dim = {'x','y','z','xy','xz','yz'};
ndim = size(scalc,2);
for d = 1:ndim
%     if sim, f1{d} = fTFEA.(dim{d}); end;            % given (FEA)
    f2{d} = scatteredInterpolant(spos,scalc(:,d));  % calculated
end
if sim
    for d=1:length(dim)
        f1{d} = fTFEA.(dim{d});
    end
end

if inp_T
    for i = 1:2
        for d=1:3
           f3{d,i} = frFEA(i).(dim{d});
        end
    end
end


eval = table();

% Evaluate stress field in the substrate (maximum principal stress)
%**************************************************************************
npos = size(spos,1);

if ndim == 6 && sim
    
    % Interpolation of the given stress field
    s_subst_given = zeros(npos,ndim);
    for d = 1:ndim
        s_subst_given(:,d) = f1{d}(spos);
    end
    % Get calculated stress field
    s_subst_calc = scalc;
    
    % Calculate maximum principal stress (MPS)
    MPS_given = zeros(npos,1);
    MPS_calc = zeros(npos,1);
    for i = 1:npos
        % Rebuild stress tensors
        sigma_given = Nye2Tensor(s_subst_given(i,:));
        sigma_calc = Nye2Tensor(s_subst_calc(i,:));
        % Calculate MPS
        MPS_given(i) = max(eig(sigma_given));
        MPS_calc(i) = max(eig(sigma_calc));
    end
    
    % Get goodness of fit (GOF) measures
    eval.substrate = {GoodnessOfFit(MPS_given,MPS_calc,spos,'MPS_')};
    fprintf('--> Substrate stress field (ROI: %.0fx%.0fx%.0f µm): MAPE = %.3f%%, MdAPE = %.3f%%, TIC = %.3f\n',...
        max(spos)-min(spos),eval.substrate{:}.MAPE,eval.substrate{:}.MdAPE,eval.substrate{:}.TIC)
    
    % Attach stress vectors to the results
    eval.substrate{1}.data{1}.sgiven = s_subst_given;
    eval.substrate{1}.data{1}.scalc = s_subst_calc;
    
elseif ndim == 6 && ~sim
    
    % Get calculated stress field
    s_subst_calc = scalc;
    
    % Calculate maximum principal stress (MPS)
    MPS_calc = zeros(size(spos,1),1);
    for i = 1:npos
        % Rebuild stress tensors
        sigma_calc = Nye2Tensor(s_subst_calc(i,:));
        % Calculate MPS
        MPS_calc(i) = max(eig(sigma_calc));
    end
    
    % Store stresses without GOF measures
    data = table(spos,MPS_calc,s_subst_calc,'VariableNames',{'pos','Tcalc','MPScalc'});
    eval.substrate = {table({data},'VariableNames',{'data'})};
    fprintf('--> Substrate stress field (DOI: %.0fx%.0fx%.0f µm)\n',max(spos)-min(spos))
    
end
        
% Evaluate tractions at cell face centroids
%**************************************************************************

% Get cell face centroids and normals
fvcell = imgparam.fvcell;
facecent = fvcell.UserData.facecent;
normals = fvcell.UserData.normals;
nfaces = size(fvcell.faces,1);

% tt = find(facecent(:,1)<8 & facecent(:,1)>-8);

if sim
    
    % Interpolation of the stress fields on cell face centroids
    s_cell_given = zeros(nfaces,ndim);
    s_cell_calc = zeros(nfaces,ndim);
    for d = 1:ndim
%         s_cell_given(:,d) = f1{d}(facecent);
        s_cell_calc(:,d) = f2{d}(facecent);
    end
    for d = 1:length(dim)
        s_cell_given(:,d) = f1{d}(facecent);
    end
    
    % Determine traction vectors at each cell face
    if ndim == 6
        T_given = zeros(nfaces,3);
        T_calc = zeros(nfaces,3);
        for i = 1:nfaces
            % Rebuild stress tensors
            sigma_given = Nye2Tensor(s_cell_given(i,:));
            sigma_calc = Nye2Tensor(s_cell_calc(i,:));
            % Calculate traction vectors (check sign of normals for method)
            T_given(i,:) = (sigma_given * -normals(i,:)')';
            T_calc(i,:) = (sigma_calc * sgn* normals(i,:)')';
            ind = {1:nfaces}';
        end
        if inp_T
            T_given_full = T_given;
            clear T_given ind
            for i = 1:2
                fa_range{i} = [min(f3{1,i}.Points(:,1)) max(f3{1,i}.Points(:,1)); ...
                    min(f3{1,i}.Points(:,2)) max(f3{1,i}.Points(:,2)); ...
                    min(f3{1,i}.Points(:,3)) max(f3{1,i}.Points(:,3))];
                
                inBounds = bsxfun(@ge, facecent, fa_range{i}(:,1)' ) & bsxfun(@le, facecent, fa_range{i}(:,2)' );
                ind{i} = find(inBounds(:,1) & inBounds(:,2) & inBounds(:,3));
                fa_facecent{i} = facecent(ind{i},:);
                
                for d = 1:3
                    input{i}(:,d) = f3{d,i}(fa_facecent{i});
                end
            end
            % Combine everything together
            T_given = vertcat(input{1},input{2});
            pos_given = vertcat(fa_facecent{:});
            T_calc_full = T_calc;
            clear T_calc
            T_calc = T_calc_full(vertcat(ind{:}),:);
            clear facecent
            facecent = pos_given;
        end
    else % Traction vectors already available
        T_given = zeros(nfaces,3);
        for i = 1:nfaces
            % Rebuild stress tensors
            sigma_given = Nye2Tensor(s_cell_given(i,:));
            % Calculate traction vectors
            T_given(i,:) = (sigma_given * -normals(i,:)')';
        end
        T_calc = s_cell_calc;
        ind = {1:nfaces}';
    end
    
    % Get goodness of fit (GOF) measures
    eval.cell = {GoodnessOfFit(T_given,T_calc,facecent,'T')};
    fprintf('--> Cell surface traction stresses (%d faces): MAPE = %.3f%%, MdAPE = %.3f%%, TIC = %.3f\n',...
        nfaces,eval.cell{:}.MAPE,eval.cell{:}.MdAPE,eval.cell{:}.TIC)
    
    % Attach stress vectors to the results
    eval.cell{1}.data{1}.sgiven = s_cell_given(vertcat(ind{:}),:);
    eval.cell{1}.data{1}.scalc = s_cell_calc(vertcat(ind{:}),:);
    
    eval.cell{1}.T_given_full = {T_given_full}; 
    eval.cell{1}.facecent_full = {fvcell.UserData.facecent}; 
    eval.cell{1}.T_calc_full = {T_calc_full};    
else
    
    % Interpolation of the calculated stress field on cell face centroids
    s_cell_calc = zeros(nfaces,ndim);
    for d = 1:ndim
        s_cell_calc(:,d) = f2{d}(facecent);
    end
    
    % Determine traction vectors at each cell face
    if ndim == 6
        T_calc = zeros(nfaces,3);
        for i = 1:nfaces
            % Rebuild stress tensors
            sigma_calc = Nye2Tensor(s_subst_calc(i,:));
            % Calculate traction vectors
            T_calc(i,:) = (sigma_calc * normals(i,:)')';
        end
    else % Traction vectors already available
        T_calc = s_cell_calc;
    end
    
    % Store stresses without GOF measures
    data = table(facecent,T_calc,'VariableNames',{'pos','Tcalc'});
    eval.substrate = {table({data},'VariableNames',{'data'})};
    fprintf('--> Cell surface traction stresses (%d faces)\n',nfaces)
    
end

end