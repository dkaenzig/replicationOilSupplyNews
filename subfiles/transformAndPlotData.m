%% Transform and plot the data

% Transform data if requested
data_diff = nan(size(data_level,1)-1,nvar);
for ii = 1:nvar
    if diffInd_diff(ii)==0
        data_diff(:,ii) = data_level(2:end,ii);
    elseif diffInd_diff(ii)==1
        data_diff(:,ii) = diff(data_level(:,ii)); %*12; % annualize
    end
end

% exogenous series
dataExo_level = dataExoRaw;
dataExo_diff  = dataExoRaw(2:end,:);

% Trim series to the relevant sample
datesString_diff = datesStringRaw(2:end);

data_level = data_level(find(strcmp(datesStringRaw,smplStart)):find(strcmp(datesStringRaw,smplEnd)),:);
data_diff  = data_diff(find(strcmp(datesString_diff,smplStart)):find(strcmp(datesString_diff,smplEnd)),:);
dataExo_level = dataExo_level(find(strcmp(datesStringRaw,smplStart)):find(strcmp(datesStringRaw,smplEnd)),:);
dataExo_diff  = dataExo_diff(find(strcmp(datesString_diff,smplStart)):find(strcmp(datesString_diff,smplEnd)),:);

sampleDates_level    = datesStringRaw(find(strcmp(datesStringRaw,smplStart)):find(strcmp(datesStringRaw,smplEnd)));
sampleDatesNum_level = datesNumRaw(find(strcmp(datesStringRaw,smplStart)):find(strcmp(datesStringRaw,smplEnd)));
sampleDates_diff     = datesString_diff(find(strcmp(datesString_diff,smplStart)):find(strcmp(datesString_diff,smplEnd)));
sampleDatesNum_diff  = datesNumRaw(find(strcmp(datesString_diff,smplStart)):find(strcmp(datesString_diff,smplEnd)));

% select relevant data
if estDiff
    data = data_diff;
    diffInd = diffInd_diff;
    dataExo = dataExo_diff;
    datesString = datesString_diff;
    sampleDates = sampleDates_diff;
    sampleDatesNum = sampleDatesNum_diff;
    varNames = varNames_level;
else
    data = data_level;
    diffInd = diffInd_level;
    dataExo = dataExo_level; 
    datesString = datesStringRaw;
    sampleDates = sampleDates_level;
    sampleDatesNum = sampleDatesNum_level;
    varNames = varNames_level;
end

% plot the data
if plotData
    figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13); 
    for ii = 1:size(data,2)
        h(ii) = subplot(2,ceil(nvar/2),ii);
        hold on
        plot(sampleDatesNum,data(:,ii),'LineWidth',1.5)
        pylims = get(gca,'ylim');
        if ~ismember(0,pylims) && ~(pylims(1)>0 || pylims(2)<0)
            l1 = line(get(gca,'xlim'),[0 0],'Color','k');%,'LineStyle','--')
            uistack(l1,'bottom');
        end
        title(varNames_paper{ii})
        %axis tight
        xlim([sampleDatesNum_level(1) sampleDatesNum_level(end)])
        grid on
        box on
    end
    pause(0.001)
    if mod(nvar,2)~=0
        pos = get(h,'Position');
        new = mean(cellfun(@(v)v(1),pos(1:2)));
        set(h(ii-1),'Position',[(pos{1}(1)+pos{2}(1))/2 pos{end}(2:end)])
        set(h(ii),'Position',[(pos{2}(1)+pos{3}(1))/2 pos{end}(2:end)])
    end
    tightfig;
    if saveFigs
        print('-dpdf', gcf, strcat(savePath,'plotData_',figName));  
    end
end