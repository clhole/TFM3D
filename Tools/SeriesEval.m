function [ serieseval ] = SeriesEval( varargin )
%SERIESEVAL summarizes the results of a series of simulations performed by
%           the function <TFM3D.m>
%Input:
%           There are two modes of using this function:
%
%           1) Add results of a dataset to the global results
%              Ex: serieseval = SeriesEval(serieseval_in,dataset)
%              <serieseval_in>  global results structure
%              <dataset>        evaluated current dataset
%              <mode>           'u' evaluate only displacement data
%                               'T' evaluate traction and displacement data
%
%           2) Summarize all results of a series after the simulation is done
%              Ex: serieseval = SeriesEval(foldername)
%              <foldername>     folder containing all datasets of the series
%              <mode>           'u' evaluate only displacement data
%                               'T' evaluate traction and displacement data
%Output:
%           <serieseval>        structure containing the summarized results  
%CL

% Parse inputs
nvarargs = length(varargin);
if nvarargs == 3
    continuous = true;                        % *** Mode 1 ***
    serieseval = varargin{1};                 % Previous results
    dataset = varargin{2};                    % Current dataset
    mode = varargin{3};                       % Data mode
    ndatasets = 1;                            % Number of datasets to add
    if isempty(serieseval)
        newstruct = true;                     % Create global results structure
        serieseval.ndatasets = 1;             % Total number of datasets
    else
        newstruct = false;                    % Don't create global results structure
        serieseval.ndatasets = serieseval.ndatasets + 1; % Total number of datasets
    end
elseif nvarargs == 2
    continuous = false;                       % *** Mode 2 ***
    serieseval = [];                          % No previous results
    foldername = varargin{1};                 % Folder containing all datasets
    mode = varargin{2};                       % Data mode
    files = dir(sprintf('.\\Results\\%s\\dataset*.mat',foldername));
    ndatasets = length(files);                % Number of datasets to add
    newstruct = true;                         % Create global results structure
    serieseval.ndatasets = ndatasets;         % Total number of datasets
end

% Repeat for all datasets to be added to the results structure
for j = 1:ndatasets
    
    % Load previous results
    if ~continuous
        dataset = load(sprintf('.\\Results\\%s\\dataset%03d.mat',foldername,j));
    end
    
    % Get parameters of the dataset
    imgpairlbl = fieldnames(dataset);
    nimgpairs = length(imgpairlbl);
    
    % Repeat for all image pairs of the current dataset
    for k = 1:nimgpairs
               
        % Get displacement calculation parameters of the current image pair
        uparam = dataset.(imgpairlbl{k}).uparam;
        ualg = uparam.ualg;
        nualg = length(ualg);
        for ua = 1:nualg
            % Get labels and number of displacement evaluation methods and parameter sets
            fieldlbl = dataset.(imgpairlbl{k}).u.(ualg{ua}).Properties.VariableNames;
            ueval{ua} = fieldlbl(length(fieldnames(uparam.(ualg{ua})))+1:end);
            nueval(ua) = size(ueval{ua},2);
            nuparamsets(ua) = size(dataset.(imgpairlbl{k}).u.(ualg{ua}),1);
            uparamlbl{ua} = fieldlbl(1:length(fieldnames(uparam.(ualg{ua}))));
        end
        
        % Get traction calculation parameters of the current image pair
        if strcmp('T',mode)
            Tparam = dataset.(imgpairlbl{k}).Tparam;
            Talg = cell(length(Tparam.Talg),1);
            for ua = 1:nualg
                for ue = 1:nueval(ua)
                    % Get labels and number of traction evaluation methods and parameter sets
                    lblTalg = fieldnames(dataset.(imgpairlbl{k}).T.(ualg{ua}).(ueval{ua}{ue}){1});
                    Talg(1:length(lblTalg),ue) = lblTalg;
                    nTalg(ua,ue) = length(lblTalg);
                    up = 1; % Consider only first parameter set to get the labels                   
                    for Ta = 1:nTalg(ua,ue)
                        fieldlbl = dataset.(imgpairlbl{k}).T.(ualg{ua}).(ueval{ua}{ue}){up}.(Talg{Ta,ue}).Properties.VariableNames;
                        Teval{Ta} = fieldlbl(length(fieldnames(Tparam.(Talg{Ta,ue})))+1:end);
                        nTeval(Ta) = size(Teval{Ta},2);
                        nTparamsets(Ta) = size(dataset.(imgpairlbl{k}).T.(ualg{ua}).(ueval{ua}{ue}){up}.(Talg{Ta,ue}),1);
                        Tparamlbl{Ta} = fieldnames(Tparam.(Talg{Ta,ue}))';
                    end
                end
            end
        end
        
        % Create new structures and tables for the results
        if newstruct
            % Create table for image generation parameters
            serieseval.imgparam = cell2table(cell(nimgpairs,4));
            serieseval.imgparam.Properties.VariableNames = {'cbeads','rbead','vxsize','PSFparam'};
            serieseval.imgparam.Properties.RowNames = imgpairlbl;
            % Create MdAPE tables
            serieseval.u.MdAPE = cell2table(cell(nueval*nuparamsets',3+nimgpairs));
            serieseval.u.MdAPE.Properties.VariableNames = {'Algorithm','ParameterSet','Evaluation',imgpairlbl{:,:}};
            if strcmp('T',mode)
                serieseval.T.MdAPE = cell2table(cell(nueval*nuparamsets'*nTeval*nTparamsets',3+nimgpairs));
                serieseval.T.MdAPE.Properties.VariableNames = {'Algorithm','ParameterSet','Evaluation',imgpairlbl{:,:}};
            end
        end
        
        % Write image generation parameters to the results file
        serieseval.imgparam.cbeads{k,1} = dataset.(imgpairlbl{k}).imgparam.cbeads;
        serieseval.imgparam.rbead{k,1} = dataset.(imgpairlbl{k}).imgparam.rbead;
        serieseval.imgparam.vxsize{k,1} = dataset.(imgpairlbl{k}).imgparam.vxsize;
        serieseval.imgparam.PSFparam{k,1} = dataset.(imgpairlbl{k}).imgparam.PSFparam;

        % Row counters of the MdAPE tables
        row_u = 0;
        row_T = 0;
        
        addrow_u = true;
        addrow_T = true;
        
        for ua = 1:nualg           
            for ue = 1:nueval(ua)
                for up = 1:nuparamsets(ua)
                    
                    if ~isfield(serieseval.u,imgpairlbl{k})
                        % Write displacement results of the very first dataset to the results structure
                        serieseval.u.(imgpairlbl{k}) = dataset.(imgpairlbl{k}).u;
                        addrow_u = false;
                    elseif addrow_u
                        % Attach displacement results of the current dataset to the results structure
                        serieseval.u.(imgpairlbl{k}).(ualg{ua}).(ueval{ua}{ue}){up} =...
                            [serieseval.u.(imgpairlbl{k}).(ualg{ua}).(ueval{ua}{ue}){up};...
                                dataset.(imgpairlbl{k}).u.(ualg{ua}).(ueval{ua}{ue}){up}];
                    end

                    % Attach displacement results to the MdAPE table
                    row_u = row_u + 1;
                    if k == 1 % 1st image pair --> create row
                        uparamval = table2cell(dataset.(imgpairlbl{k}).u.(ualg{ua})(up,1:size(uparamlbl{ua},2)));
                        uparamstr = [uparamlbl{ua};cellfun(@(x) strrep(['[' strrep(num2str(x,'%g,'),' ',''),' '],', ',']'),uparamval,'Unif',false)];
                        serieseval.u.MdAPE.Algorithm{row_u} = ualg{ua};
                        serieseval.u.MdAPE.ParameterSet{row_u} = [sprintf('%d)',up),' ',sprintf('%s=%s ',uparamstr{:})];
                        serieseval.u.MdAPE.Evaluation{row_u} = ueval{ua}{ue};
                    end
                    serieseval.u.MdAPE.(imgpairlbl{k}){row_u} = nanmean(serieseval.u.(imgpairlbl{k}).(ualg{ua}).(ueval{ua}{ue}){up}.MdAPE);

                    if strcmp('T',mode)
                        for Ta = 1:nTalg(ua,ue)
                            
                            for Te = 1:nTeval(Ta)
                                for Tp = 1:nTparamsets(Ta)

                                    if ~isfield(serieseval.T,imgpairlbl{k})
                                        % Write traction results of the very first dataset to the results structure
                                        serieseval.T.(imgpairlbl{k}) = dataset.(imgpairlbl{k}).T;
                                        addrow_T = false;
                                    elseif addrow_T
                                        % Attach traction results of the current dataset to the results structure
                                         serieseval.T.(imgpairlbl{k}).(ualg{ua}).(ueval{ua}{ue}){up}.(Talg{Ta,ue}).(Teval{Ta}{Te}){Tp} =...
                                            [serieseval.T.(imgpairlbl{k}).(ualg{ua}).(ueval{ua}{ue}){up}.(Talg{Ta,ue}).(Teval{Ta}{Te}){Tp};...
                                                dataset.(imgpairlbl{k}).T.(ualg{ua}).(ueval{ua}{ue}){up}.(Talg{Ta,ue}).(Teval{Ta}{Te}){Tp}];
                                    end

                                    % Attach traction results to the MdAPE table
                                    row_T = row_T + 1;
                                    if k == 1 % 1st image pair --> create row
                                        serieseval.T.MdAPE.Algorithm{row_T} = sprintf('%s, %s',ualg{ua},Talg{Ta,ue});
                                        Tparamval = table2cell(dataset.(imgpairlbl{k}).T.(ualg{ua}).(ueval{ua}{ue}){up}.(Talg{Ta,ue})(Tp,1:size(Tparamlbl{Ta},2)));
                                        Tparamstr = [Tparamlbl{Ta};cellfun(@(x) strrep(['[' strrep(num2str(x,'%g,'),' ',''),' '],', ',']'),Tparamval,'Unif',false)];
                                        serieseval.T.MdAPE.ParameterSet{row_T} = [sprintf('%d%s)',up,char('a'+Tp-1)),' ',sprintf('%s=%s ',uparamstr{:}),sprintf('%s=%s ',Tparamstr{:})];
                                        serieseval.T.MdAPE.Evaluation{row_T} = sprintf('%s, %s',ueval{ua}{ue},Teval{Ta}{Te});
                                    end
                                    try 
                                        serieseval.T.MdAPE.(imgpairlbl{k}){row_T} = nanmean(serieseval.T.(imgpairlbl{k}).(ualg{ua}).(ueval{ua}{ue}){up}.(Talg{Ta,ue}).(Teval{Ta}{Te}){Tp}.MdAPE);
                                    catch
                                        serieseval.T.MdAPE.(imgpairlbl{k}){row_T} = NaN;
                                    end

                                end %Tp
                            end %Te
                        end %Ta
                    end

                end %up
            end %ue
        end %ua
        
        newstruct = false;
        
    end %k
    
end %j

% Delete empty rows in MdAPE tables
if strcmp('T',mode)
    ind = find(~cellfun(@isempty,serieseval.T.MdAPE{:,1}));
    serieseval.T.MdAPE = serieseval.T.MdAPE(ind,:);
end

end

