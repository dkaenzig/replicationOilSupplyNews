%% Replication files for "The macroeconomic effects of oil supply news"
% This is the main file to produce the main results in the paper

% All code was written and tested in Matlab R2019b
% The results are saved in ../results of the current directory
% Auxiliary functions (tools) are located in /auxfiles
% Subfunctions are located in /subfiles

% Diego R. Känzig
% LBS, December 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% add tools and subfiles directories
addpath(genpath('auxfiles'))
addpath(genpath('subfiles'))

% initialize random number generator to be able to replicate results exactly
rng default

% Set text interpreter for figures to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Read in data

load data/OilDataM

Tall = size(datesNumRaw,1);
dataFrequency = 'M'; 

% create deterministic variables
const = ones(Tall,1);
trend1 = (1:Tall)';
trend2 = trend1.^2;
zlbDummy = datesNumRaw>=datesNumRaw(strcmp(datesStringRaw,'2008M09'));
gmDummy = datesNumRaw>=datesNumRaw(strcmp(datesStringRaw,'1984M01'));
seasDummies = zeros(Tall,11);
for id = 1:11
    seasDummies(id:12:end,id) = 1;
end

%% Settings
% settings for VAR model

% select the estimation sample
smplStart = '1974M01'; 
smplEnd   = '2017M12'; 

% select the identification sample
% range has to be contained in estimation sample
smplStartProxy = '1975M01';  
smplEndProxy   = '2017M12'; 
% we censor missing values to zero and use the same sample, adjusted for
% lags. Alternatively, we can use '1983M04', producing very similar results

% normalize trends
trend1 = trend1./trend1(strcmp(datesStringRaw,smplStart));
trend2 = trend2./trend2(strcmp(datesStringRaw,smplStart));

% VAR specifics
p          = 12;          % lag order
dataExoRaw = [const]; 
estDiff    = false;       % VAR specification (estimate in differences or levels)
estType    = 'proxy';     % Estimation: external instrument VAR ('proxy') 
                          % or heteroskedasticity-based ('hetero')
horizon    = 50;          % horizon for IRFs
shockType  = 'custom';    % one standard deviation 'sd' or 'custom'
shockSize  = 10;          % if custom, specify shock size here
alpha      = 0.1;         % Significance level for bands (alpha=0.1 => 90% CIs (two SD))
alpha2     = 0.32;
nsim       = 1000;        % number of simulations in bootstrap (in paper 10000 are used)
bootType   = 'mbb1block'; % Options: mbb bootstrap 'mbb1block', iid bootstrap 'iid',
                          % wild bootstrap 'wild' (not valid but added for comparison)

% proxy
ncontract = 15;       % select futures contract 
% options: front (1) to 12-month contract (13), spot (14), PC(1-12) (15), PC(spot+0-12) (16)

% switches
plotData   = false;  % plot the data
doProxyVAR = true;   % run proxy VAR
doLagCrit  = false;  % compute lag length criteria
getFEVDs   = false;  % get FEVDs
getHD      = true;   % get HD
doLP       = true;   % also compute IRFs using LP on shock

saveFigs   = true;       % save figures to disk

%% Data manipulations
% select data
data_level = [log(POIL)*100-log(CPI/100)*100 log(OILPROD)*100 log(OILSTOCKS)*100 log(WORLDIP)*100 log(IP)*100 log(CPI)*100];  
% set series names
varNames_level = {'POIL','OILPROD','OILSTOCKS','WORLDIP','IP','CPI'}; 
varNames_diff  = varNames_level;
varNames_paper = {'Real oil price','World oil production','World oil inventories','World industrial production','U.S. industrial production','U.S. CPI'};
varNames_paperVD = {'Real oil price','Oil production','Oil inventories','World IP','U.S. IP','U.S. CPI'};

% number of variables in VAR
nvar = size(data_level,2);    

% set integration order (only relevant if estimated in differences)
diffInd_level = zeros(1,nvar);
diffInd_diff  = [0 1 1 1 1 1];

% select relevant sample, transform the data if requested and plot the data
transformAndPlotData;

%% Load instrument

loadProxy;

% plot proxy
figure('DefaultAxesFontSize',13); 
hold on
sampleDatesNumSel = sampleDatesNum(smplStartProxyVARInd:smplEndProxyVARInd);
plot(sampleDatesNumSel,proxy)
plot(sampleDatesNumSel(abs(proxy)>7 & sampleDatesNumSel<2015),proxy(abs(proxy)>7 & sampleDatesNumSel<2015),'ro')
xlim([sampleDatesNum(smplStartProxyVARInd) sampleDatesNum(smplEndProxyVARInd)])
ylim([-15 15])
xlim([1984 2017])
text(1987,10.5,'5 Aug 1986','fontsize',11)
text(1996,-8.6,'14 Nov 2001','fontsize',11)
text(2009,-11.25,'27 Nov 2014','fontsize',11)
ylabel('Revision in oil price expectations [\%]')
line(get(gca,'xlim'),[0 0],'Color','k')
grid on
box on
if saveFigs
    saveas(gcf,'Figures/proxyPlot','epsc2')
end


%% Estimate model: either external instrument/proxy VAR or heteroskedasticity-based VAR

if strcmp(estType,"proxy")
    runProxyVAR; 
elseif strcmp(estType,"hetero")
    runRigobonVAR;
end

%% Plot impulse responses

time = (0:horizon)';    % time horizon for IRFs and FEVDs

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
for j=1:nvar 
    h(j) = subplot(ceil(nvar/3),3,j);

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_pe(1,j); IRFslower_pe(1:end,j); flipud([IRFsupper_pe(1:end,j); IRFslower_pe(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 

    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_pe(1,j); IRFslower2_pe(1:end,j); flipud([IRFsupper2_pe(1:end,j); IRFslower2_pe(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    p1=plot(time, IRFs_pe(:,j),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    box on
    grid on ;hold off;
    title(varNames_paper{j}) 
    ylabel('\%');
    xlim([0,horizon]);
    if dataFrequency == 'M'
        xlabel('Months');
        xticks(0:10:horizon)
    elseif dataFrequency == 'Q'
        xlabel('Quarters');
        xticks(0:5:horizon)
    end   
end
pause(0.001)
h=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
set(h,'Visible','off');
if k==1 && strcmp(estType,'proxy')
    string_1stage = ['First stage regression: F: ',num2str(olsEst.r1.F,' %2.2f'),', robust F: ',num2str(olsEst.r1.Frobust,' %2.2f'),', $R^2$: ',num2str(olsEst.r1.R2*100,' %1.2f'),'\%, Adjusted $R^2$: ',num2str(olsEst.r1.R2adj*100,' %1.2f'),'\%'];
    text('Position',[-0.16 -0.002],'string',string_1stage,'FontSize',14);
end
tightfig;
if saveFigs
    print('-dpdf', gcf, strcat('figures/IRFs_',estType));  
end


%% Historical decomposition of oil price

if getHD
    % Historical decomposition
    varCompan = varxcompan(varEst);    % companion form
    F = varCompan.F;                   % max(abs(eig(F))): VAR is stable 
    V       = oilSupplyNewsShock';     % structural errors 
    nvarXeq = nvar*p;                  % number of lagged endogenous per equation

    % Compute historical decompositions

    % Contribution of the shock
    A0_big = zeros(nvarXeq,nvar);
    if strcmp(shockType,'custom')
        A0_big(1:nvar,1) = b1unit;
    else
        A0_big(1:nvar,1) = b1;
    end
    Icomp = [eye(nvar) zeros(nvar,(p-1)*nvar)];
    HDshock_big = zeros(p*nvar,T+1);
    HDshock_temp = zeros(nvar,T+1);
    V_big = zeros(nvar,T+1); % matrix of shocks conformable with companion
    V_big(1,2:end) = V;
    for i = 2:T+1
        HDshock_big(:,i) = A0_big*V_big(:,i) + F*HDshock_big(:,i-1);
        HDshock_temp(:,i) =  Icomp*HDshock_big(:,i);
    end
    HDshock = HDshock_temp(:,2:end)';
    
    % bootstrap it 
    bootHDshock = zeros(T,nvar,nsim);
    for isim = 1:nsim
        varBoot = struct;
        varBoot.B = bootBs(:,:,isim);
        varBoot.Sigma = Sigma;
        varBoot.nexo = nexo;
        varBootCompan = varxcompan(varEst);    % companion form
        bootF = varBootCompan.F;
        bootV = bootShocks(:,isim)';                         % structural errors 
        
        % Contribution of the shock
        bootA0_big = zeros(nvarXeq,nvar);
        bootA0_big(1:nvar,1) = bootb1s(:,isim);
        Icomp = [eye(nvar) zeros(nvar,(p-1)*nvar)];
        bootHDshock_big = zeros(p*nvar,T+1);
        bootHDshock_temp = zeros(nvar,T+1);
        bootV_big = zeros(nvar,T+1); % matrix of shocks conformable with companion
        bootV_big(1,2:end) = bootV;
        for i = 2:T+1
            bootHDshock_big(:,i) = bootA0_big*bootV_big(:,i) + bootF*bootHDshock_big(:,i-1);
            bootHDshock_temp(:,i) =  Icomp*bootHDshock_big(:,i);
        end
        bootHDshock(:,:,isim) = bootHDshock_temp(:,2:end)';
    end
    
    HDshockmed = quantile(bootHDshock, 0.5, 3);
    HDshockupper = quantile(bootHDshock, 1-alpha/2, 3)-HDshockmed+HDshock;
    HDshocklower = quantile(bootHDshock, alpha/2, 3)-HDshockmed+HDshock;
    HDshockupper2 = quantile(bootHDshock, 1-alpha2/2, 3)-HDshockmed+HDshock;
    HDshocklower2 = quantile(bootHDshock, alpha2/2, 3)-HDshockmed+HDshock;
    
    timeHD = sampleDatesNum(smplStartProxyVARInd:smplEndProxyVARInd);
    oilDates = [1978+8/12 1980+9/12 1985+11/12 1990+7/12 1997+6/12 2002+10/12 2008+8/12 2014+10/12];
    
    figure('Position',[100 100 1000 350],'DefaultAxesFontSize',13)
    hold on
    hh=fill([timeHD(1); timeHD(1:end); flipud([timeHD(1:end); timeHD(end)])],[HDshockupper(1,1); HDshocklower(1:end,1); flipud([HDshockupper(1:end,1); HDshocklower(end,1)])]-mean(HDshock(:,1)),[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 
    
    hh=fill([timeHD(1); timeHD(1:end); flipud([timeHD(1:end); timeHD(end)])],[HDshockupper2(1,1); HDshocklower2(1:end,1); flipud([HDshockupper2(1:end,1); HDshocklower2(end,1)])]-mean(HDshock(:,1)),[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none'); 
    set(hh,'edgealpha',.4);
    f1=plot(sampleDatesNum(smplStartProxyVARInd:smplEndProxyVARInd),varEst.Y(:,1)-mean(varEst.Y(:,1)),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2,'LineStyle',':');
    f2=plot(sampleDatesNum(smplStartProxyVARInd:smplEndProxyVARInd),HDshock(:,1)-mean(HDshock(:,1)),'Color','k','LineWidth',1.5);
    for ii=1:length(oilDates)
       xline(oilDates(ii),'LineWidth',1.2); 
    end
    grid on
    box on
    xlim([1975 2017])
    ylim([-150 150])
    ylabel('\%')
    tightfig;
    
    fs = 10;
    legend([f1 f2],'Real oil price','Contribution oil supply news','Location','Southeast','Interpreter','latex')
    str={'Iranian', 'revolution'};
    annotation('textbox',[0.077, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
    str={'Iran-Iraq','war'};
    annotation('textbox',[0.2, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
    str={'OPEC','collapse'};
    annotation('textbox',[0.31, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
    str={'Gulf war'};
    annotation('textbox',[0.414, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
    str={'Asian fin.','crisis'};
    annotation('textbox',[0.49, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
    str={'Venezuelan','crisis'};
    annotation('textbox',[0.596, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
    str={'Global fin.','crisis'};
    annotation('textbox',[0.73, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
    str={'Oil crash','2014'};
    annotation('textbox',[0.872, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);

    if saveFigs
        print('-dpdf', gcf, strcat('figures/HD_',estType));
    end
    
end


%% Compute IRFs using local projections to shock

pLP = 1;
horizonLP = horizon;

keepSampleFixed = true;

colval = [0.8500, 0.3250, 0.0980]; 

% run LPs on shock
if doLP

    IRFs_LP = zeros(horizonLP+1,nvar);

    for hh = 0:horizonLP
        for ii = 1:nvar
            if keepSampleFixed
                yi = data(smplStartProxyVARInd+hh:end-horizonLP+hh,ii);

                Xr = oilSupplyNewsShock(1:end-horizonLP);  
                for jj = 1:pLP
                    Xr = [Xr data(smplStartProxyVARInd-jj:end-horizonLP-jj,ii)];
                end
                Xr = [Xr dataExo(smplStartProxyInd:smplEndProxyInd-horizonLP,2:end)];
            else
                yi = data(smplStartProxyVARInd+hh:end,ii);

                Xr = oilSupplyNewsShock(1:end-hh);  
                for jj = 1:pLP
                    Xr = [Xr data(smplStartProxyVARInd-jj:end-hh-jj,ii)];
                end
                Xr = [Xr dataExo(smplStartProxyInd:smplEndProxyInd-hh,2:end)];
            end
            
            olsLP = olsest(Xr,yi,true);
            IRFs_LP(hh+1,ii) = olsLP.bhat(1);
        end
    end
    IRFs_LP = IRFs_LP./IRFs_LP(1,1)*shockSize;
    
    % bands
    IRFsboot_LP = zeros(horizonLP+1,nvar,nsim);
    
    for j = 1:nsim
        for hh = 0:horizonLP
            for ii = 1:nvar
                if keepSampleFixed
                    yi = bootDatas(smplStartProxyVARInd+hh:end-horizonLP+hh,ii,j);

                    Xr = bootShocks(1:end-horizonLP,j);  
                    for jj = 1:pLP
                        Xr = [Xr bootDatas(smplStartProxyVARInd-jj:end-horizonLP-jj,ii,j)];
                    end
                    Xr = [Xr dataExo(smplStartProxyInd:smplEndProxyInd-horizonLP,2:end)];
                else
                    yi = bootDatas(smplStartProxyVARInd+hh:end,ii,j);

                    Xr = bootShocks(1:end-hh,j);  
                    for jj = 1:pLP
                        Xr = [Xr bootDatas(smplStartProxyVARInd-jj:end-hh-jj,ii,j)];
                    end
                    Xr = [Xr dataExo(smplStartProxyInd:smplEndProxyInd-hh,2:end)];
                end

                olsLP = olsest(Xr,yi);
                IRFsboot_LP(hh+1,ii,j) = olsLP.bhat(1);
            end
        end
        IRFsboot_LP(:,:,j) = IRFsboot_LP(:,:,j)./IRFsboot_LP(1,1,j)*IRFs_LP(1,1);
    end
    
    % rescale bootstrapped IRFs (center around sample estimates) and get quantiles
    IRFsmed_LP = quantile(IRFsboot_LP, 0.5, 3);

    IRFsupper_LP = quantile(IRFsboot_LP, 1-alpha/2, 3)-IRFsmed_LP+IRFs_LP;  % rescaling does not affect the ordering, thus it is fine to take quantile first
    IRFslower_LP = quantile(IRFsboot_LP, alpha/2, 3)-IRFsmed_LP+IRFs_LP;

    IRFsupper2_LP = quantile(IRFsboot_LP, 1-alpha2/2, 3)-IRFsmed_LP+IRFs_LP;  % rescaling does not affect the ordering, thus it is fine to take quantile first
    IRFslower2_LP = quantile(IRFsboot_LP, alpha2/2, 3)-IRFsmed_LP+IRFs_LP;
    
    time = (0:horizonLP)';

    % Plot it
    figure('Position',[10 10 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
    if IRFs_LP(1,1)>0
        signIRFs = 1;
    else
        signIRFs = -1;
    end
    for j=1:nvar %variable
        subplot(2,ceil(nvar/2),j) %note the index
        
        hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_LP(1,j); IRFslower_LP(1:end,j); flipud([IRFsupper_LP(1:end,j); IRFslower_LP(end,j)])],[0.1, 0.4470, 0.7410]); 
        set(hh,'facealpha',.2);
        set(hh,'edgecolor','none'); 

        hold on;
        hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_LP(1,j); IRFslower2_LP(1:end,j); flipud([IRFsupper2_LP(1:end,j); IRFslower2_LP(end,j)])],[0.1, 0.4470, 0.7410]); 
        set(hh,'facealpha',.4);
        set(hh,'edgecolor','none');
        
        p2=plot(time, IRFs_pe(1:horizonLP+1,j), 'Color',colval, 'Linewidth', 1.5,'LineStyle','-');
        p1=plot(time, IRFs_LP(:,j),'k', 'Linewidth', 1.5); 
        plot(time, IRFsupper_pe(1:horizonLP+1,j), 'Color',colval,'Linestyle','--');
        plot(time, IRFslower_pe(1:horizonLP+1,j), 'Color',colval,'Linestyle','--');
        plot(time, IRFsupper2_pe(1:horizonLP+1,j), 'Color',colval,'Linestyle',':');
        plot(time, IRFslower2_pe(1:horizonLP+1,j), 'Color',colval,'Linestyle',':');
        grid on ;hold off;
        if j==1
            legend('LP-IV','Proxy-VAR','AutoUpdate','off')
        end
        title(varNames_paper{j}) 
        if dataFrequency == 'M'
            xlabel('Months');
        elseif dataFrequency == 'Q'
            xlabel('Quarters');
        end  
        ylabel('\%');
        if ~ismember(0,get(gca,'ylim'))
            line(get(gca,'xlim'),[0 0],'Color','k')%,'LineStyle','--')
        end
        xlim([0,horizonLP]);
        xticks([0:10:horizonLP]);
        if j==1
            legend([p1 p2],{'LP','VAR'})
        end
    end
    tightfig;
    if saveFigs
        print('-dpdf', gcf, strcat('figures/IRFs_LP_',estType));  
    end
end


