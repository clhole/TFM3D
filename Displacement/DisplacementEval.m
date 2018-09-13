function [ eval ] = DisplacementEval( upos, ucalc, ualg, imgparam, uparam )
%DISPLACEMENTEVAL evaluates the displacement field at the specified locations
%Input
%   <upos>        positions where displacements were calculated [µm]
%   <ucalc>       calculated displacements [µm]
%   <ualg>        current displacement calculation algorithm
%   <imgparam>    structure containing the following fields:
%                 .range        nodal range [µm]                           [nx3]
%                 .fvcell       cell model patch object [µm]                
%                 .fuFEA        FEA displacement interpolant [µm]
%   <uparam>      structure containing the following fields:
%                 .itp          displacement interpolation locations       {1xn}
%                               (options: 'bead','grid','cell')
%                 .itp_spacing  grid spacing [µm] for itp = 'grid'         [1xn]
%                 .ROIsize      region of interest [% cell size]           [1xn]
%Output
%   <eval>        table containing the evaluated displacement fields
%CL

% Get image parameters
range = imgparam.range;
fvcell = imgparam.fvcell;
if isfield(imgparam,'fuFEA') % For simulated data
    sim = true;
    fuFEA = imgparam.fuFEA;
else % For experimental data
    sim = false;
end

% Get displacement calculation parameters
itp = uparam.itp;
ROIsize = uparam.ROIsize;

% Determine the cell-centered region of interest (ROI) for analysis
cellrange = transpose([min(fvcell.vertices);max(fvcell.vertices)]);
nROI = size(ROIsize,2);
for i = 1:nROI
    if isnan(ROIsize(i))
        ROI{i} = range;
        ROIname{i} = 'full';
    else
        ROI{i} = cellrange*ROIsize(i)/100;
        ROIname{i} = sprintf('ROI%d',ROIsize(i));
    end
end

% Prepare displacement interpolants
dim = ['x','y','z'];
for d = 1:3
    if sim, f1{d} = fuFEA.(dim(d)); end;            % given (FEA)
    f2{d} = scatteredInterpolant(upos,ucalc(:,d));  % calculated
end

eval = table();

% Interpolate displacements on detected neutral bead positions
%**************************************************************************
if any(strcmp('bead',itp)) && strcmp('PDVC',ualg)

    % Get bead positions
    beadpos = upos;
    
    if sim
        
        % Interpolation
        for d = 1:3
            u_bead_given(:,d) = f1{d}(beadpos);  % given displacements
            u_bead_calc(:,d) = f2{d}(beadpos);   % calculated displacements
        end
        
        for i = 1:nROI
            method = sprintf('bead_%s',ROIname{i});
            % Get displacements within the region of interest (ROI)
            u_bead = [beadpos,u_bead_given,u_bead_calc];
            u_bead = u_bead(u_bead(:,1)>=ROI{i}(1,1) & u_bead(:,1)<=ROI{i}(1,2) &...
                        u_bead(:,2)>=ROI{i}(2,1) & u_bead(:,2)<=ROI{i}(2,2) &...
                        u_bead(:,3)>=ROI{i}(3,1) & u_bead(:,3)<=ROI{i}(3,2),:);
            % Calculate goodness of fit (GOF) measures
            eval.(method) = {GoodnessOfFit(u_bead(:,4:6),u_bead(:,7:9),u_bead(:,1:3),'u')};
            eval.(method){:}.ROI = ROI(i);
            fprintf('--> Bead displacement (ROI: %.0fx%.0fx%.0f µm, %d beads): MAPE = %.3f%%, TIC = %.3f\n',...
                ROI{i}(:,2)-ROI{i}(:,1),size(u_bead,1),eval.(method){:}.MAPE,eval.(method){:}.TIC)
        end
        
    else
        
        % Interpolation
        for d = 1:3
            u_bead_calc(:,d) = f2{d}(beadpos);
        end

        for i = 1:nROI
            method = sprintf('bead_%s',ROIname{i});
            % Get displacements within the region of interest (ROI)
            u_bead = [beadpos,u_bead_calc];
            u_bead = u_bead(u_bead(:,1)>=ROI{i}(1,1) & u_bead(:,1)<=ROI{i}(1,2) &...
                        u_bead(:,2)>=ROI{i}(2,1) & u_bead(:,2)<=ROI{i}(2,2) &...
                        u_bead(:,3)>=ROI{i}(3,1) & u_bead(:,3)<=ROI{i}(3,2),:);
            % Store displacements without GOF measures
            data = table(u_bead(:,1:3),u_bead(:,4:6),'VariableNames',{'pos','ucalc'});
            eval.(method) = {table({data},ROI(i),'VariableNames',{'data','ROI'})};
            fprintf('--> Bead displacement (ROI: %.0fx%.0fx%.0f µm, %d beads)\n',...
                ROI{i}(:,2)-ROI{i}(:,1),size(u_bead,1))
        end 
        
    end
    
end

% Interpolate displacements on grid nodes
%**************************************************************************
if any(strcmp('grid',itp)) || strcmp('GDVC',ualg) || strcmp('FIDVC',ualg) || strcmp('known',ualg)
    
    itp_spacing = uparam.itp_spacing;
    nitp = size(itp_spacing,2);
    
    % Repeat for each specified grid spacing
    for sp = 1:nitp
    
        % Prepare grid within the range of displacement data
        xgv = ceil(range(1,1)):itp_spacing(sp):floor(range(1,2));
        ygv = ceil(range(2,1)):itp_spacing(sp):floor(range(2,2));
        zgv = ceil(range(3,1)):itp_spacing(sp):floor(range(3,2));
        [x,y,z] = meshgrid(xgv,ygv,zgv);
        nodes = [x(:),y(:),z(:)];

        if sim

            % Interpolation
            for d = 1:3
                u_grid_given(:,d) = f1{d}(nodes);  % given displacements
                u_grid_calc(:,d) = f2{d}(nodes);   % calculated displacements
            end

            for i = 1:nROI
                method = sprintf('grid_%04.0f_%s',10*itp_spacing(sp),ROIname{i});
                % Get displacements within the region of interest (ROI)
                u_grid = [nodes,u_grid_given,u_grid_calc];
                u_grid = u_grid(nodes(:,1)>=ROI{i}(1,1) & nodes(:,1)<=ROI{i}(1,2) &...
                            nodes(:,2)>=ROI{i}(2,1) & nodes(:,2)<=ROI{i}(2,2) &...
                            nodes(:,3)>=ROI{i}(3,1) & nodes(:,3)<=ROI{i}(3,2),:);
                % Avoid evaluating data points within the cell
                inshp = inpolyhedron(fvcell,u_grid(:,1:3));
                [row,~] = find(inshp);
                u_grid(row,:) = NaN;
                % Calculate goodness of fit (GOF) measures 
                eval.(method) = {GoodnessOfFit(u_grid(:,4:6),u_grid(:,7:9),u_grid(:,1:3),'u')};
                eval.(method){:}.ROI = ROI(i);
                fprintf('--> Grid displacement (ROI: %.0fx%.0fx%.0f µm, %d nodes, spacing = %.1f µm): MAPE = %.3f%%, TIC = %.3f\n',...
                     ROI{i}(:,2)-ROI{i}(:,1),size(u_grid,1),itp_spacing(sp),eval.(method){:}.MAPE,eval.(method){:}.TIC)
            end

            clearvars u_grid_given u_grid_calc
            
        else
            
            % Interpolation
            for d = 1:3
                u_grid_calc(:,d) = f2{d}(nodes);
            end

            for i = 1:nROI
                method = sprintf('grid_%03.0f_%s',itp_spacing(sp),ROIname{i});
                % Get displacements within the region of interest (ROI)
                u_grid = [nodes,u_grid_calc];
                u_grid = u_grid(nodes(:,1)>=ROI{i}(1,1) & nodes(:,1)<=ROI{i}(1,2) &...
                            nodes(:,2)>=ROI{i}(2,1) & nodes(:,2)<=ROI{i}(2,2) &...
                            nodes(:,3)>=ROI{i}(3,1) & nodes(:,3)<=ROI{i}(3,2),:);
                % Avoid evaluating data points within the cell
                inshp = inpolyhedron(fvcell,u_grid(:,1:3));
                [row,~] = find(inshp);
                u_grid(row,:) = NaN;
                % Store displacements without GOF measures
                data = table(u_grid(:,1:3),u_grid(:,4:6),'VariableNames',{'pos','ucalc'});
                eval.(method) = {table({data},ROI(i),'VariableNames',{'data','ROI'})};
                fprintf('--> Grid displacement (ROI: %.0fx%.0fx%.0f µm, %d nodes, spacing = %.1d µm)\n',...
                    ROI{i}(:,2)-ROI{i}(:,1),size(u_grid,1),itp_spacing(sp))
            end

            clearvars u_grid_given u_grid_calc
            
        end
        
    end
    
end

% Interpolate displacements on cell vertices
%**************************************************************************
if any(strcmp('cell',itp))
    
    if sim
    
        % Interpolation
        for d = 1:3
            u_cell_given(:,d) = f1{d}(fvcell.vertices);  % given displacements
            u_cell_calc(:,d) = f2{d}(fvcell.vertices);   % calculated displacements
        end
        
        % Calculate goodness of fit (GOF) measures
        eval.cell = {GoodnessOfFit(u_cell_given,u_cell_calc,fvcell.vertices,'u')};
        eval.cell{:}.ROI = {'cell vertices'};
        fprintf('--> Cell surface displacement (%d vertices): MAPE = %.3f%%, TIC = %.3f\n',...
            size(fvcell.vertices,1),eval.cell{:}.MAPE,eval.cell{:}.TIC)
        
    else
        
        % Interpolation
        for d = 1:3
            u_cell_calc(:,d) = f2{d}(fvcell.vertices);
        end

        % Store displacements without GOF measures
        data = table(fvcell.vertices,u_cell_calc,'VariableNames',{'pos','ucalc'});
        eval.(method) = {table({data},{'cell vertices'},'VariableNames',{'data','ROI'})};
        fprintf('--> Cell surface displacement (%d vertices)\n',size(fvcell.vertices,1))
        
    end
    
end

end