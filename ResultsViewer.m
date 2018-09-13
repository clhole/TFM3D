function ResultsViewer( varargin )
%EVALTABLE       displays the results of a simulation series
%Input:
%  <foldername>  (optional) folder containing the results of a simulation series
%  <mode>        (optional)
%                'T'    display traction results (default)
%                'u'    display displacement results
%CL

% Get foldername and mode
if ~isempty(varargin)
    foldername = varargin{1};
    path = sprintf('.\\Results\\%s',foldername);
    mode = varargin{2};
else
    path = uigetdir(strcat(pwd,'\Results'),'Select Folder');
    [~,foldername,~] = fileparts(path);
    choice = questdlg('Select results to display','Mode','Displacements','Traction','Traction');
    switch choice
        case 'Displacements', mode = 'u';
        case 'Traction', mode = 'T';
    end

end

% Load series GOF data
serieseval = load(strcat(path,'\serieseval.mat'));
nimgpairs = size(serieseval.imgparam,1);
for k = 1:nimgpairs
    lbl{k,:} = sprintf('%.3f|%.3f|[%.1f,%.1f,%.1f]|[%s]',cell2mat(serieseval.imgparam{k,1:3}),...
        sprintf('%.2f,%.2f,%.2f',cell2mat(serieseval.imgparam{k,4})));
end

% Create table
%**************************************************************************
f = figure('name',foldername,'NumberTitle','off','menubar','none','visible','off');
t = uitable('Parent',f);
try
    t.Data = serieseval.(mode).MdAPE{:,:};
catch
    t.Data = serieseval.(mode).MAPE{:,:};
end
% Make columns sortable
jscrollpane = findjobj(t);
jtable = jscrollpane.getViewport.getView;
pause(0.05); % Time buffer to prevent java thread exceptions
jtable.setSortable(true);
jtable.setAutoResort(true);
jtable.setMultiColumnSortable(true);
jtable.setPreserveSelectionsAfterSorting(true);
% Set table properties
if strcmp('u',mode)
    title('Displacement Field: Averaged MdAPE of the Series');
else
    title('Traction Field: Averaged MdAPE of the Series');
end
axis off
t.FontSize = 10;
t.ColumnName = {' | |Algorithm',' | |ParameterSet',...
    ['                                                    cbead|',...
    '                                                    rbead|',...
    '                                                    vxsize|',...
    'Evaluation                        PSF[NA,ni,ns]'],lbl{:}};
t.ColumnWidth = {150,520,270,'auto'};
t.RowName = [];
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);
f.Position(3) = t.Extent(3)+40;
f.Position(4) = t.Extent(4)+80;
% Show the figure at screen center
movegui(f,'center')
set(get(0,'children'),'visible','on')

% Set up context menu for selecting plots
%**************************************************************************
c = uicontextmenu(f);
t.UIContextMenu = c;

% Create menu items for the uicontextmenu
m1 = uimenu(c,'Label','Plot bead distribution','Callback',{@PlotSelected,t,foldername,serieseval,mode});
if strcmp('T',mode)
    m2 = uimenu(c,'Label','Plot traction field','Callback',{@PlotSelected,t,foldername,serieseval,mode});
else
    m2 = uimenu(c,'Label','Plot displacement field','Callback',{@PlotSelected,t,foldername,serieseval,mode});
end
m3 = uimenu(c,'Label','Plot error distribution','Callback',{@PlotSelected,t,foldername,serieseval,mode});
m4 = uimenu(c,'Label','Bland-Altman plot','Callback',{@PlotSelected,t,foldername,serieseval,mode});
m5 = uimenu(c,'Label','Plot APE histogram','Callback',{@PlotSelected,t,foldername,serieseval,mode});
m6 = uimenu(c,'Label','Plot MAPE histogram','Callback',{@PlotSelected,t,foldername,serieseval,mode});

end

% Callback function to display selected plot
%**************************************************************************
function PlotSelected( source, event, t, foldername, serieseval, mode )

% Get selection
jscrollpane = findjobj(t);
jtable = jscrollpane.getViewport.getView;
row = jtable.getSelectedRow();
col = jtable.getSelectedColumn();

if col > 2
    
    MAPE = jtable.getModel.getValueAt(row,col).doubleValueReal;
    
    if ~isnan(MAPE)
    
        % Get data
        imgpair = sprintf('imgpair%03d',col-2);
        [ualg,remain1] = strtok(jtable.getModel.getValueAt(row,0),',');
        uparamset = sscanf(jtable.getModel.getValueAt(row,1),'%d');
        [ueval,remain2] = strtok(jtable.getModel.getValueAt(row,2),',');
        if strcmp('T',mode)
            Talg = remain1(3:end);
            alphabet = char('a'+(1:26)-1);
            c = char(sscanf(jtable.getModel.getValueAt(row,1),'%*d%c%*s'));
            Tparamset = strfind(alphabet,c);
            Teval = remain2(3:end);
        end

        % Find representative dataset for plotting (MAPE = approx. series MAPE)
        ndatasets = serieseval.ndatasets;
        if strcmp('T',mode)
            for j = 1:ndatasets
                datasetMAPE(j) = nanmean(serieseval.T.(imgpair).(ualg).(ueval){uparamset}.(Talg).(Teval){Tparamset}.data{j}.APE);
            end
        else
            for j = 1:ndatasets
                datasetMAPE(j) = nanmean(serieseval.u.(imgpair).(ualg).(ueval){uparamset}.data{j}.APE);
            end
        end
        [~,datasetnumber] = min(abs(datasetMAPE-MAPE));

        % Display selected plot
        hbox1 = msgbox('Loading data...','modal');
        switch source.Label
            case 'Plot bead distribution'
                dataset = load(sprintf('.\\Results\\%s\\dataset%03d.mat',foldername,datasetnumber),imgpair);
                ROI = serieseval.u.(imgpair).(ualg).(ueval){uparamset,1}.ROI{j,1};
                PlotBeadImagePair(dataset.(imgpair),ROI,datasetnumber) 
            case 'Plot displacement field'
                dataset = load(sprintf('.\\Results\\%s\\dataset%03d.mat',foldername,datasetnumber),imgpair);
                PlotDisplacement(dataset.(imgpair),ualg,uparamset,ueval,datasetnumber);
            case 'Plot traction field'
                dataset = load(sprintf('.\\Results\\%s\\dataset%03d.mat',foldername,datasetnumber),imgpair);
                PlotTraction(dataset.(imgpair),ualg,uparamset,ueval,Talg,Tparamset,Teval,datasetnumber);
            case 'Plot error distribution'
                if strcmp('T',mode)
                    PlotErrorDistribution(serieseval,imgpair,ualg,uparamset,ueval,Talg,Tparamset,Teval);
                else
                    PlotErrorDistribution(serieseval,imgpair,ualg,uparamset,ueval);
                end
            case 'Bland-Altman plot'
                if strcmp('T',mode)
                    BlandAltmanPlot(serieseval,imgpair,ualg,uparamset,ueval,Talg,Tparamset,Teval);
                else
                    BlandAltmanPlot(serieseval,imgpair,ualg,uparamset,ueval);
                end
            case 'Plot APE histogram'
                if strcmp('T',mode)
                    PlotAPEHistogram(serieseval,imgpair,ualg,uparamset,ueval,Talg,Tparamset,Teval);
                else
                    PlotAPEHistogram(serieseval,imgpair,ualg,uparamset,ueval);
                end
            case 'Plot MAPE histogram'
                if strcmp('T',mode)
                    PlotMAPEHistogram(serieseval,imgpair,ualg,uparamset,ueval,Talg,Tparamset,Teval);
                else
                    PlotMAPEHistogram(serieseval,imgpair,ualg,uparamset,ueval);
                end
        end
        close(hbox1);
        
    end
    
else
    
    hbox2 = msgbox('Please select a MAPE value.','modal');
    
end

end

