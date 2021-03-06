function PlotErrorDistribution( serieseval, imgpairlbl, ualg ,uparamset, ueval, varargin )
%PlotErrorDistribution displays the error distribution of a simulation
%                      series depending on location and magnitude
%Input:
%   <imgpaireval>      image pair series evaluation generated by <TFMSim3D.m>
%   <imgpairlbl>       image pair name
%   <ualg>             name of displacement calculation algorithm
%   <uparamset>        parameter set number for displacement calculation
%   <ueval>            name of displacement evaluation method
%   <Talg>             (optional) name of traction calculation algorithm
%   <Tparamset>        (optional) parameter set number for traction calculation
%   <Teval>            (optional) name of traction evaluation method
%CL

% Determine method
if ~isempty(varargin)
    % Use traction data
    method = 'T';
    Talg = varargin{1};
    Tparamset = varargin{2};
    Teval = varargin{3};
    if strcmp('substrate',Teval)
        Tdata = 'MPS_given';
    	lbl = 'MPS';
    else
        Tdata = 'Tgiven';
        lbl = 'T';
    end
    unit = 'Pa';
    title_T = sprintf(', T: %s (%s)',Talg,strrep(Teval,'_',', '));
else
    % Use displacement data
    method = 'u';
    lbl = 'u';
    unit = '�m';
    title_T = ' ';
end

ndatasets = serieseval.ndatasets;
c = lines(ndatasets); % Colors

figure
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.7 0.8])
h = suptitle({'Absolute percentage error (APE) distribution',...
    sprintf('u: %s (%s)%s',ualg,strrep(ueval,'_',', '),title_T),' '});
set(h,'FontSize',11,'FontWeight','bold')

for j = 1:ndatasets
    
    if strcmp('T',method)
        loc = serieseval.T.(imgpairlbl).(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.data{j}.pos;
        mag = serieseval.T.(imgpairlbl).(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.data{j}.(Tdata);
        APE = serieseval.T.(imgpairlbl).(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.data{j}.APE;
    else
        loc = serieseval.u.(imgpairlbl).(ualg).(ueval){uparamset}.data{j}.pos;
        mag = serieseval.u.(imgpairlbl).(ualg).(ueval){uparamset}.data{j}.ugiven;
        APE = serieseval.u.(imgpairlbl).(ualg).(ueval){uparamset}.data{j}.APE;
    end
    
    subplot(2,2,1)
    hold on
    scatter(loc(:,1),APE,'x','MarkerEdgeColor',c(j,:)); % APE vs x
    xlabel('x [�m]')
    ylabel('APE [%]')
    
    subplot(2,2,2)
    hold on
    scatter(loc(:,2),APE,'x','MarkerEdgeColor',c(j,:)); % APE vs y
    xlabel('y [�m]')
    ylabel('APE [%]')
    lstr{j} = sprintf('imgpair%03d',j);
    legend(lstr);
    
    subplot(2,2,3)
    hold on
    scatter(loc(:,3),APE,'x','MarkerEdgeColor',c(j,:)); % APE vs z
    xlabel('z [�m]')
    ylabel('APE [%]')
    
    subplot(2,2,4)
    hold on
    scatter(round(1000*sqrt(sum(mag.^2,2)))/1000,APE,'x','MarkerEdgeColor',c(j,:)); % APE vs magnitude
    xlabel(sprintf('known |%s| [%s]',lbl,unit))
    ylabel('APE [%]')
    
end

end

