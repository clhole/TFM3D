function [ GOF ] = GoodnessOfFit( given, calc, varargin )
%GOODNESSOFFIT determines goodness of fit parameters
%Input:
%   <given>    given values
%   <calc>     calculated values
%   <pos>      (optional) position where values are compared
%   <mode>     (optional) mode for labeling, e.g. 'u' or 'T'
%Output:
%   <GOF>      table with GOF results
%CL

% Number of data points
npoints = size(given,1);

% Parse optional inputs
nvarargs = length(varargin);
if nvarargs == 2
    pos = varargin{1};
    mode = varargin{2};
elseif nvarargs == 1
    if ischar(varargin{1})
        pos = NaN(npoints,3);
        mode = varargin{1};
    else
        pos = varargin{1};
        mode = '';
    end
else
    pos = NaN(npoints,3);
    mode = '';
end

% Preallocation
DM = zeros(npoints,1);
DA = zeros(npoints,1);
EVM = zeros(npoints,1);
APE = zeros(npoints,1);

for i = 1:npoints
    % Deviation of magnitude (DM)
    DM(i) = norm(given(i,:)) - norm(calc(i,:));
    % Deviation of angle (DA)
    DA(i) = acos((given(i,:)*calc(i,:)')/(norm(given(i,:))*norm(calc(i,:))))*180/pi;
    % Error vector magnitude (EVM)
    EVM(i) = norm(given(i,:)-calc(i,:));
    % Modified absolute percentage error (APE)
    APE(i) = 100*EVM(i)/max(norm(given(i,:)),norm(calc(i,:)));
    % Absolute deviation of magnitude (ADM)
    ADM(i) = DM(i)/norm(given(i,:));
end

% Create data table
data = [table(pos),table(given,calc,APE,DM,DA,...
    'VariableNames',{[mode,'given'];[mode,'calc'];'APE';'DM';'DA'})];

% Mean, Median and standard deviation of DM and DA
meanDM = nanmean(DM);
stdDM = nanstd(DM);
medianDM = nanmedian(DM);
meanDA = nanmean(DA);
stdDA = nanstd(DA);
medianDA = nanmedian(DA);

% Mean absolute percentage error (MAPE)
MAPE = nanmean(APE);
% Median absolute percentage error (MAPE)
MdAPE = nanmedian(APE);
% Theil's inequality coefficient (TIC)
TIC = sqrt(nanmean(EVM.^2))/(sqrt(nanmean(sum(given.^2,2)))+sqrt(nanmean(sum(calc.^2,2))));

% Create GOF table
GOF = table({data},meanDM,stdDM,meanDA,stdDA,MAPE,MdAPE,TIC,...
    'VariableNames',{'data','meanDM','stdDM','meanDA','stdDA','MAPE','MdAPE','TIC'});
    
end
