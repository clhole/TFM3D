function PlotDisplacement( imgpair, ualg, uparamset, uevalmethod, varargin )
%PLOTDISPLACEMENT displays the calculated displacement fields
%Input:
%  <imgpair>        structure array containing the fields:
%                   imgparam.fvcell          cell model patch object [µm]
%                   u.(method).bead.pos0     detected neutral bead positions [µm]
%                   u.(method).bead.calc     bead-based displacement field [µm]
%                   u.(method).grid.fcalc    grid-based displacement interpolant [µm]
%                   u.(method).cell.calc     cell vertex-based displacement field [µm]
%  <ualg>           name of displacement calculation algorithm
%  <uparamset>      parameter set number for displacement calculation
%  <evalmethod>     name of displacement evaluation method
%  <datasetnumber>  (optional) datasetnumber to display in the figure title
%CL

% If handed over, display the number of the dataset in the figure title    
if ~isempty(varargin)
    datasetnumber = varargin{1};
    reprdataset = sprintf('Representative dataset: No. %03d',datasetnumber);
else
    reprdataset = '';
end

% Get parameters
fvcell = imgpair.imgparam.fvcell;
itp = strtok(uevalmethod,'_');

% Bead-based displacement field
if strcmp('bead',itp)

    % Get data
    pos = imgpair.u.(ualg).(uevalmethod){uparamset}.data{1}.pos;
    ucalc = imgpair.u.(ualg).(uevalmethod){uparamset}.data{1}.ucalc;
    ROI = imgpair.u.(ualg).(uevalmethod){uparamset}.ROI{1};
    MAPE = imgpair.u.(ualg).(uevalmethod){uparamset}.MAPE;
    TIC = imgpair.u.(ualg).(uevalmethod){uparamset}.TIC;
    
    figure
    hold on
    % Plot displacement field
    nbeads = size(pos,1);
    if nbeads < 1000 % Plot colored displacement vectors (slower)
        c = jet(64);
        ncolors = length(c);
        uabs = sqrt(sum(ucalc.^2,2));
        umax = max(uabs);
        umin = min(uabs);
        for n = 1:size(pos,1)
            cind = ceil((uabs(n)-umin)/(umax-umin)*(ncolors-1)+1);
            quiver3(pos(n,1),pos(n,2),pos(n,3),ucalc(n,1),ucalc(n,2),ucalc(n,3),2,...
                'Color',[c(cind,1),c(cind,2),c(cind,3)],'LineWidth',0.8,'MaxHeadSize',0.7)
        end
        colormap(jet)
        hcb = colorbar;
        caxis([umin,umax]);
        xlabel(hcb,'Displacement [µm]','FontSize',12);
    else % Black displacement vectors (faster)
        quiver3(pos(:,1),pos(:,2),pos(:,3),ucalc(:,1),ucalc(:,2),ucalc(:,3),5,'k')
    end
    % Plot cell
    patch(fvcell,'FaceColor','g','FaceAlpha',0.3)
    % Set figure properties
    set(gcf,'units','normalized','outerposition',[0.12 0.08 0.72 0.78])
    title({sprintf('%s',reprdataset);...
        ['Displacement field interpolated on detected neutral bead positions; '...
        sprintf('DOI: %.1fx%.1fx%.1f µm, %d beads',ROI(:,2)-ROI(:,1),size(pos,1))];...
        sprintf('Method: %s',ualg);...
        sprintf('MAPE = %.3f%%, TIC = %.3f',MAPE,TIC)})
    xlabel('x [µm]')
    ylabel('y [µm]')
    zlabel('z [µm]')
    axis equal
    axis([ROI(1,:),ROI(2,:),ROI(3,:)]);
    grid on
    view(3)

end

% Grid-based displacement field
if strcmp('grid',itp)
    
    % Get data
    pos = imgpair.u.(ualg).(uevalmethod){uparamset}.data{1}.pos;
    ucalc = imgpair.u.(ualg).(uevalmethod){uparamset}.data{1}.ucalc;
    ROI = imgpair.u.(ualg).(uevalmethod){uparamset}.ROI{1};
    MAPE = imgpair.u.(ualg).(uevalmethod){uparamset}.MAPE;
    TIC = imgpair.u.(ualg).(uevalmethod){uparamset}.TIC;
    
    % Downsampling for clearer vector plot
    sp = 50; % Spacing between vectors to plot
    pos_s = pos(1:sp:end,:);
    ucalc_s = ucalc(1:sp:end,:);
    
    % Prepare displacement data for volumetric slice plot
    [x,y,z] = ReconstructMeshgrid(pos); % Data points
    uabs = reshape(sqrt(sum(ucalc.^2,2)),size(x)); % Magnitudes
    uabs(isnan(uabs)) = 0;
    nslices = 10;
    zslice = linspace(min(pos(:,3)),max(pos(:,3)),nslices);

    figure
    hold on
    % Plot traction field
    quiver3(pos_s(:,1),pos_s(:,2),pos_s(:,3),ucalc_s(:,1),ucalc_s(:,2),ucalc_s(:,3),1,'k')
    h = slice(x,y,z,uabs,[0],[0],[0]); % Volumetric slice plot
    set(h, 'EdgeColor','none', 'FaceColor','interp')
    alpha(.2)
    colormap(jet)
    hcb = colorbar;
    % Plot cell
    patch(fvcell,'FaceColor','g','FaceAlpha',0.3)
    % Set figure properties
    set(gcf,'units','normalized','outerposition',[0.14 0.06 0.74 0.76])
%     title({sprintf('%s',reprdataset);...
%         ['Displacement field interpolated on a regular grid; '...
%         sprintf('DOI: %.1fx%.1fx%.1f µm, %d nodes',ROI(:,2)-ROI(:,1),size(pos,1))];...
%         sprintf('Method: %s',ualg);...
%         sprintf('MAPE = %.3f%%, TIC = %.3f',MAPE,TIC)})
    xlabel(hcb,'Displacement [µm]','FontSize', 20);
    set(gca,'XTick', [-40 0 40],'FontSize',20)
    set(gca,'YTick', [-40 0 40],'FontSize',20)
    set(gca,'ZTick', [-40 0 40],'FontSize',20)
    xlabel('x [µm]')
    ylabel('y [µm]')
    zlabel('z [µm]')
    axis equal
    axis([ROI(1,:),ROI(2,:),ROI(3,:)]);
    grid on
    view(3)
    print(strcat(pwd,'\figures\displacement_grid'),'-dtiff','-r0')

end

% Cell vertex-based displacement field
if strcmp('cell',itp)

    % Get data
    pos = imgpair.u.(ualg).(uevalmethod){uparamset}.data{1,1}.pos;
    ucalc = imgpair.u.(ualg).(uevalmethod){uparamset}.data{1,1}.ucalc;
    uabs = sqrt(sum(ucalc.^2,2)); % Magnitude
    MAPE = imgpair.u.(ualg).(uevalmethod){uparamset}.MAPE;
    TIC = imgpair.u.(ualg).(uevalmethod){uparamset}.TIC;
    
    % Downsampling for clearer vector plot
    sp = 2; % Spacing between vectors to plot
    vert = pos(1:sp:end,:);
    ucalc_s = ucalc(1:sp:end,:);

    figure
    hold on
    % Plot displacement field
    quiver3(vert(:,1),vert(:,2),vert(:,3),ucalc_s(:,1),ucalc_s(:,2),ucalc_s(:,3),0,'k')
    % Plot cell
    patch(fvcell,'FaceVertexCData',uabs,'FaceColor','interp','FaceAlpha',0.8,'EdgeColor','none')
    colormap(jet)
    hcb = colorbar;
    % Set figure properties
    set(gcf,'units','normalized','outerposition',[0.16 0.04 0.76 0.74])
    title({sprintf('%s',reprdataset);...
        'Displacement field interpolated on cell surface vertices';...
        sprintf('Method: %s',ualg);...
        sprintf('MAPE = %.3f%%, TIC = %.3f',MAPE,TIC)})
    xlabel(hcb,'Displacement [µm]','FontSize', 12);
    xlabel('x [µm]')
    ylabel('y [µm]')
    zlabel('z [µm]')
    axis equal
    grid on
    view(3)

end

end
    

