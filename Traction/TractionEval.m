function [ eval ] = TractionEval( spos, scalc, imgparam )
%TRACTIONEVAL   evaluates the stress field in the substrate (maximum
%               principal stress) and on the cell surface ('tractions')
%Input:
%   <spos>      positions where the stress field was calculated [�m]
%   <scalc>     stress tensors in Nye notation [sx,sy,sz,sxy,sxz,syz]
%               OR traction vectors [Tx,Ty,Tz]
%   <imgparam>  structure containing the following fields:
%               .fvcell>  cell model patch object incl. face centroids [�m]
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

% Prepare stress field interpolants
dim = {'x','y','z','xy','xz','yz'};
ndim = size(scalc,2);
for d = 1:ndim
    if sim, f1{d} = fTFEA.(dim{d}); end;            % given (FEA)
    f2{d} = scatteredInterpolant(spos,scalc(:,d));  % calculated
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
    fprintf('--> Substrate stress field (ROI: %.0fx%.0fx%.0f �m): MAPE = %.3f%%, TIC = %.3f\n',...
        max(spos)-min(spos),eval.substrate{:}.MAPE,eval.substrate{:}.TIC)
    
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
    fprintf('--> Substrate stress field (DOI: %.0fx%.0fx%.0f �m)\n',max(spos)-min(spos))
    
end
        
% Evaluate tractions at cell face centroids
%**************************************************************************

% Get cell face centroids and normals
fvcell = imgparam.fvcell;
facecent = fvcell.UserData.facecent;
normals = fvcell.UserData.normals;
nfaces = size(fvcell.faces,1);

tt = find(facecent(:,1)<8 & facecent(:,1)>-8);

if sim
    
    % Interpolation of the stress fields on cell face centroids
    s_cell_given = zeros(nfaces,ndim);
    s_cell_calc = zeros(nfaces,ndim);
    for d = 1:ndim
        s_cell_given(:,d) = f1{d}(facecent);
        s_cell_calc(:,d) = f2{d}(facecent);
    end
    
    % Determine traction vectors at each cell face
    if ndim == 6
        T_given = zeros(nfaces,3);
        T_calc = zeros(nfaces,3);
        for i = 1:nfaces
            % Rebuild stress tensors
            sigma_given = Nye2Tensor(s_cell_given(i,:));
            sigma_calc = Nye2Tensor(s_cell_calc(i,:));
            % Calculate traction vectors
            T_given(i,:) = (sigma_given * -normals(i,:)')';
            T_calc(i,:) = (sigma_calc * -normals(i,:)')';
        end
    else % Traction vectors already available 
        T_given = s_cell_given;
        T_calc = s_cell_calc;
    end
%     T_given(tt,:) = [];
%     T_calc(tt,:) = [];
%     facecent(tt,:) = [];
    
    % Get goodness of fit (GOF) measures
    eval.cell = {GoodnessOfFit(T_given,T_calc,facecent,'T')};
    fprintf('--> Cell surface traction stresses (%d faces): MAPE = %.3f%%, TIC = %.3f\n',...
        nfaces,eval.cell{:}.MAPE,eval.cell{:}.TIC)
    
    % Attach stress vectors to the results
    eval.cell{1}.data{1}.sgiven = s_cell_given;
    eval.cell{1}.data{1}.scalc = s_cell_calc;
    
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