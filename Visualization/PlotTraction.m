function PlotTraction( imgpair, ualg, uparamset, ueval, Talg, Tparamset, Teval, varargin )
%PLOTTRACTION plots the calculated traction stresses
%Input:
%  <imgpair>        structure array containing the fields:
%                   imgparam.fvcell             cell model patch object [µm]
%                   T.(method).substrate.Tcalc  substrate stress field [Pa]
%                   T.(method).cell.Tcalc       cell surface stress [Pa]
%  <ualg>           name of displacement calculation algorithm
%  <uparamset>      parameter set number for displacement calculation
%  <uevalmethod>    name of displacement evaluation method
%  <Talg>           name of traction calculation algorithm
%  <Tparamset>      parameter set number for traction calculation
%  <Tevalmethod>    name of traction evaluation method
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

% Substrate maximum principal stress (MPS)
if strcmp('substrate',Teval)

    % Get data
    pos = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.data{1}.pos;
    MPS_calc = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.data{1}.MPS_calc;
    MAPE = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.MAPE;
    TIC = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.TIC;
    
    % Prepare traction data for volumetric slice plot
    spacing = 1;
    range = [min(pos);max(pos)]';
    xgv = ceil(range(1,1)*10)/10:spacing:floor(range(1,2)*10)/10;
    ygv = ceil(range(2,1)*10)/10:spacing:floor(range(2,2)*10)/10;
    zgv = ceil(range(3,1)*10)/10:spacing:floor(range(3,2)*10)/10;
    [x,y,z] = meshgrid(xgv,ygv,zgv);
    % Interpolate MPS on grid nodes
    f = scatteredInterpolant(pos,MPS_calc);
    MPS = f(x,y,z);
    nslices = 10;
    zslice = linspace(min(pos(:,3)),max(pos(:,3)),nslices);

    figure
    hold on
    % Plot MPS field
    h = slice(x,y,z,MPS,[],[],zslice); % Volumetric slice plot
    set(h, 'EdgeColor','none', 'FaceColor','interp')
    alpha(.2)
    colormap(jet)
    hcb = colorbar;
    % Plot cell
    patch(fvcell,'FaceColor','g','FaceAlpha',0.3)
    % Set figure properties
    set(gcf,'units','normalized','outerposition',[0.14 0.06 0.74 0.76])
    title({sprintf('%s',reprdataset);...
        ['Maximum principal stress (MPS); '...
        sprintf('%d nodes',size(pos,1))];...
        sprintf('Methods (u,T): %s, %s',ualg,Talg);...
        sprintf('MAPE = %.3f%%, TIC = %.3f',MAPE,TIC)})
    xlabel(hcb,'MPS [Pa]','FontSize', 12);
    xlabel('x [µm]')
    ylabel('y [µm]')
    zlabel('z [µm]')
    axis equal
    grid on
    view(3)

end

% Tractions on cell vertices
if strcmp('cell',Teval)

    % Get data
    pos = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.data{1}.pos;
    Tcalc = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.data{1}.Tcalc;
%     pos_full = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.facecent_full{1};
%     Tcalc_full = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.T_calc_full{1};
    for d = 1:3
        fTFEA = scatteredInterpolant(pos,Tcalc(:,d));
        Tcalc_itp(:,d) = fTFEA(fvcell.vertices); % Interpolate tractions to vertices for plotting
    end
    Tgiven = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.data{1}.Tgiven;
    Tabs = sqrt(sum(Tcalc_itp.^2,2)); % Magnitude
    MAPE = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.MAPE;
    TIC = imgpair.T.(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.TIC;

    f = figure;
    hold on
    % Plot traction field
    quiver3(pos(:,1),pos(:,2),pos(:,3),Tcalc(:,1),Tcalc(:,2),Tcalc(:,3),3,'k')
%     quiver3(pos(:,1),pos(:,2),pos(:,3),Tgiven(:,1),Tgiven(:,2),Tgiven(:,3),3,'r')
    % Plot cell
    colormap(jet);
    patch(fvcell,'FaceVertexCData',Tabs,'FaceColor','interp','FaceAlpha',1)
    hcb = colorbar;
    % Set figure properties
    set(gcf,'units','normalized','outerposition',[0.16 0.04 0.76 0.74])
%     title({sprintf('%s',reprdataset);...
%         'Tractions on cell vertices';...
%         sprintf('Methods (u,T): %s, %s',ualg,Talg);...
%         sprintf('MAPE = %.3f%%, TIC = %.3f',MAPE,TIC)})
    xlabel(hcb,'Traction [Pa]','FontSize', 18);
    axis equal
    axis([-15 15 -15 15 -15 15]);
    set(gca,'XTick', [-10 0 10],'FontSize',20)
    set(gca,'YTick', [-10 0 10],'FontSize',20)
    set(gca,'ZTick', [-10 0 10],'FontSize',20)
    xlabel('x [µm]')
    ylabel('y [µm]')
    zlabel('z [µm]')
    grid on
    view(3)
    print(f,strcat(pwd,'\figures\traction_cell'),'-dtiff','-r0')

end

end
    

