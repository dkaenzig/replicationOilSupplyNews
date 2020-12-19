%% Run heteroskedasticity-based VAR a la Rigobon/Nakamura-Steinsson

% load placebo
load data/OilSurprisesMLogControl

proxyControlRaw = [oilPlaceboWTIM(:,ncontract)]; 

proxyControl = proxyControlRaw(smplStartProxyInd:smplEndProxyInd,:);   % we loose the first p values in the VAR

% variance ratio
vratio = var(oilProxiesWTI(:,ncontract))/var(oilPlaceboWTI(:,ncontract));


% Brown-Forsythe test
temp = [[oilProxiesWTI(:,ncontract); nan(length(oilPlaceboWTI(:,ncontract))-length(oilProxiesWTI(:,ncontract)),1)] oilPlaceboWTI(:,ncontract)];
[pvalue,Fstat] = vartestn(temp,'TestType','BrownForsythe','Display','off'); % 'LeveneAbsolute'); %


% plot proxy and placebo
smplStartProxyIndPlot = find(strcmp(sampleDatesProxy,'1984M04'));
smplEndProxyIndPlot   = find(strcmp(sampleDatesProxy,smplEndProxy));
smplStartProxyVARIndPlot = find(strcmp(sampleDates,'1984M04'));
smplEndProxyVARIndPlot   = find(strcmp(sampleDates,smplEndProxy));
    
% compute pdfs
[fproxy,xproxy,bw] = ksdensity(oilProxiesWTI(:,ncontract));
x = xproxy;
pd_kernel = fitdist(oilProxiesWTI(:,ncontract),'kernel','Kernel','epanechnikov'); % 
feproxy = pdf(pd_kernel,x);

[fcontrol,xcontrol,bw] = ksdensity(oilPlaceboWTI(:,ncontract));
pd_kernel = fitdist(oilPlaceboWTI(:,ncontract),'kernel','Kernel','epanechnikov','Bandwidth',0.45); % 
fecontrol = pdf(pd_kernel,xcontrol);

% figure
figure('Position',[100 100 900 350],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
subplot(1,2,1)
hold on
sampleDatesNumSel = sampleDatesNum(smplStartProxyVARIndPlot:smplEndProxyVARIndPlot);
p2=plot(sampleDatesNumSel,proxyControlRaw(smplStartProxyIndPlot:smplEndProxyIndPlot),'color',[0.8500, 0.3250, 0.0980]);
p1=plot(sampleDatesNumSel,proxyRaw(smplStartProxyIndPlot:smplEndProxyIndPlot),'color',[0, 0.4470, 0.7410]);
xlim([sampleDatesNum(smplStartProxyVARIndPlot) sampleDatesNum(smplEndProxyVARIndPlot)])
ylim([-15 15])
xticks([1985:5:2017])
ylabel('\%')
line(get(gca,'xlim'),[0 0],'Color','k')
title('Level')
grid on
box on
%legend([p1 p2],{'Proxy','Control'})
str = {strcat('$V_{ann.}/V_{control}: ',sprintf('%2.2f',vratio),'$'), ...
       'H$_0$: $V_{ann.} \leq V_{control}$', ...
       strcat('F-stat: $',sprintf('%2.2f',Fstat.fstat),'$')};
annotation('textbox',[0.14, 0.23, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',11);
subplot(1,2,2)
hold on
h1=plot(xproxy,feproxy,'LineWidth',1.5);
h2=plot(xcontrol,fecontrol,'LineWidth',1.5);
grid on 
box on
title('PDF')
xlabel('$x$')
ylabel('$f(x)$')
legend([h1 h2],'Announcement','Control','Autoupdate','off','Interpreter','latex')
uistack(h1,'top')
ylim([0 0.45])
xlim([-15 15])
if saveFigs
    saveas(gcf,'figures/proxyPlaceboPlot','epsc2')
end


% run reduced-form VAR 
varEst = varxest(data,dataExo,p);

% identification using heteroskedasticity in the data
nexo = size(dataExo,2);
U = varEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation  

% Split data into treatment and control sample
futuresStart = find(strcmp(sampleDatesProxy(smplStartProxyInd:smplEndProxyInd), '1983M04'));
statementMindSel = statementMind(smplStartProxyInd:smplEndProxyInd,1);
statementMindPlaceboSel = placeboMind(smplStartProxyInd:smplEndProxyInd,1);
indsR1 = logical(statementMindSel);
indsR2 = logical(statementMindPlaceboSel);

% IV-implementaion
T_OPEC = size(proxy(indsR1,:),1);
T_Control = size(proxyControl(indsR2,:),1);

XrIV = [(proxy(indsR1,:)-mean(proxy(indsR1,:)))/sqrt(T_OPEC); (proxyControl(indsR2,:)-mean(proxyControl(indsR2,:)))/sqrt(T_Control)];
ZrIV = [(proxy(indsR1,:)-mean(proxy(indsR1,:)))/sqrt(T_OPEC); -(proxyControl(indsR2,:)-mean(proxyControl(indsR2,:)))/sqrt(T_Control)];
yiIV = [(U(indsR1,:)-mean(U(indsR1,:)))/sqrt(T_OPEC); (U(indsR2,:)-mean(U(indsR2,:)))/sqrt(T_Control)];

% first stage
olsEst = olsest(ZrIV,XrIV,true,true);
uhat = olsEst.yhat;
            
% second stage
b21ib11_2SLS    =   [uhat]\yiIV;  
b1 = b21ib11_2SLS';      % 2 SLS coefficients               
b1unit = b1./b1(1,1);

% compute IRFs to shock
IRFs_pe = varirfsingle(varEst.B(:,1+nexo:end),b1,p,horizon,diffInd);
if strcmp(shockType,'custom')
    IRFs_pe = IRFs_pe./IRFs_pe(1,1)*shockSize;
end


% compute the confidence bands using bootstrapping
bootIRFs = nan(horizon+1,nvar,nsim);
bootb1s = nan(nvar,nsim);
bootBs = nan(nvar,size(varEst.B,2),nsim);
bootShocks = nan(T,nsim);
bootDatas = zeros(T+p,nvar,nsim); 

ProxyCount = zeros(nsim,size(proxy,2));
T_est = varEst.T; % length of estimation sample

if strcmp(bootType,'mbb1block')
    % if identification sample is shorter that estimation sample, censor
    % unobserved values to zero
    proxyLong = zeros(T_est, size(proxy,2));
    proxyLong(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:) =  proxy;

    proxyControlLong = zeros(T_est, size(proxyControl,2));
    proxyControlLong(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:) =  proxyControl;
    
    indsR1Long = zeros(T_est, size(indsR1,2));
    indsR1Long(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:) =  indsR1;
    indsR2Long = zeros(T_est, size(indsR2,2));
    indsR2Long(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:) =  indsR2;
    
    BlockSize = round(5.03*T_est^0.25);
    nBlock = ceil(T_est/BlockSize);
    VARBlocks = zeros(BlockSize,nvar,T_est-BlockSize+1);
    ProxyBlocks = zeros(BlockSize,size(proxyLong,2),T_est-BlockSize+1);
    ProxyControlBlocks = zeros(BlockSize,size(proxyLong,2),T_est-BlockSize+1);
    indsR1Blocks = zeros(BlockSize,size(indsR1Long,2),T_est-BlockSize+1);
    indsR2Blocks = zeros(BlockSize,size(indsR2Long,2),T_est-BlockSize+1);
    for j = 1:T_est-BlockSize+1
        VARBlocks(:,:,j) = varEst.U(j:BlockSize+j-1,:);
        ProxyBlocks(:,:,j) = proxyLong(j:BlockSize+j-1,:);
        ProxyControlBlocks(:,:,j) = proxyControlLong(j:BlockSize+j-1,:);
        indsR1Blocks(:,:,j) = indsR1Long(j:BlockSize+j-1,:);
        indsR2Blocks(:,:,j) = indsR2Long(j:BlockSize+j-1,:);
    end

    % center the bootstrapped VAR errors
    VARcentering = zeros(BlockSize,nvar);
    for j = 1:BlockSize
        VARcentering(j,:) = mean(varEst.U(j:T_est-BlockSize+j,:),1);
    end
    VARcentering = repmat(VARcentering,[nBlock,1]);
    VARcentering = VARcentering(1:T_est,:);

    %center the bootstrapped proxy variables
    Proxycentering = zeros(BlockSize,size(proxyLong,2));
    ProxyControlcentering = zeros(BlockSize,size(proxyControlLong,2));
    for j = 1:BlockSize 
        subProxy = proxyLong(j:T_est-BlockSize+j,:);
        subProxyControl = proxyControlLong(j:T_est-BlockSize+j,:);
        %Proxycentering(j,:) = mean(subProxy((subProxy(:,1) ~= 0),1),1); 
        % account for non-zero mean instrument:
        Proxycentering(j,:) = mean(subProxy((subProxy(:,1) ~= 0),1),1) - mean(proxyLong((proxyLong(:,1) ~= 0),1),1);
        ProxyControlcentering(j,:) = mean(subProxyControl((subProxyControl(:,1) ~= 0),1),1) - mean(proxyControlLong((proxyControlLong(:,1) ~= 0),1),1);

    end
    Proxycentering = repmat(Proxycentering,[nBlock,1]);
    Proxycentering = Proxycentering(1:T_est,:);
    
    ProxyControlcentering = repmat(ProxyControlcentering,[nBlock,1]);
    ProxyControlcentering = ProxyControlcentering(1:T_est,:);
    
end

bootindsR1 = indsR1;
bootindsR2 = indsR2;
j = 1;
while j <= nsim
    % generate artificial data
    
    if strcmp(bootType,'wild')
        % Wild bootstrap using Rademacher (should not use because invalid)
        
        rr = randn(T_est,1); %1-2*(rand(T_est,1)>0.5); %      % use wild bootstrap based on Rademacher distribution
        bootU = (varEst.U).*(rr*ones(1,nvar)); % draw from reduced-form shocks
        bootProxy = proxy.*(rr(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:)*ones(1,size(proxy,2))); % draw from proxy (note the adjustment because of shorter sample)
        bootProxyControl = proxyControl.*(rr(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:)*ones(1,size(proxy,2))); % draw from proxy (note the adjustment because of shorter sample)

        % simulate VAR starting from initial values
        Xexo = varEst.Xexo;
        
        bootData = zeros(T_est+p,nvar); 
        bootData(1:p,:) = data(1:p,:); %initial values of y, same for all j
        for i = p+1:T_est+p
            bootData(i,:)= varEst.B*[Xexo(i-p,:)'; vec(fliplr(bootData(i-p:i-1,:)'))] ...
                             + bootU(i-p,:)'; % bootstrap
        end
        
        % re-estimate the VAR
        bootvarEst = varxest(bootData,dataExo,p);
        
        % only use proxy sample for identification
        bootU = bootvarEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation  
        bootSigma = bootU'*bootU/(T-p*nvar-nexo);
        
    elseif strcmp(bootType,'iid')
        
        % Because of iid assumption, we can proceed in two steps: 1. first
        % draw with replacement from full VAR residual sample to generate
        % artificial data, 2. draw with replacement from bootstrapped VAR residuals 
        % and the proxy/placebo
        
        % 1. Draw from full VAR residual sample
        index = randsample(1:T_est,T_est,true)';    % test code using (1:T)';
        bootU = varEst.U(index,:);
        
        % simulate VAR starting from initial values
        Xexo = varEst.Xexo;
        
        bootData = zeros(T_est+p,nvar); 
        bootData(1:p,:) = data(1:p,:); %initial values of y, same for all j
        for i = p+1:T_est+p
            bootData(i,:)= varEst.B*[Xexo(i-p,:)'; vec(fliplr(bootData(i-p:i-1,:)'))] ...
                             + bootU(i-p,:)'; % bootstrap
        end
        
        % re-estimate the VAR
        bootvarEst = varxest(bootData,dataExo,p);
        
        % get residuals for identification sample
        bootUhat = bootvarEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation
        
        % 2. Draw from bootstrapped residuals and proxy on identification
        % sample
        index2 = randsample(1:T,T,true)';
        bootProxy = proxy(index2,:);
        bootProxyControl = proxyControl(index2,:);
        bootU = bootUhat(index2,:);
        
        % count the number proxy variables not censored to zero
        ProxyCount(j,:) = sum(abs(bootProxy) > 0,1);
        
        if ProxyCount(j,:)<15
            continue
        end  
        
    elseif strcmp(bootType,'mbb1block')
        % Moving block bootstrap (Lundsford and Jentsch) using one block
        % type
        
        %draw bootstrapped residuals and proxies
        index = ceil((T_est - BlockSize + 1)*rand(nBlock,1));
        bootU = zeros(nBlock*BlockSize,nvar);
        for kk = 1:nBlock
            bootU(1+BlockSize*(kk-1):BlockSize*kk,:) = VARBlocks(:,:,index(kk,1));
        end
        bootU = bootU(1:T_est,:);
        
        bootProxy = zeros(nBlock*BlockSize,size(proxy,2));
        bootProxyControl = zeros(nBlock*BlockSize,size(proxyControl,2));
        for kk = 1:nBlock
            bootProxy(1+BlockSize*(kk-1):BlockSize*kk,:) = ProxyBlocks(:,:,index(kk,1));
            bootProxyControl(1+BlockSize*(kk-1):BlockSize*kk,:) = ProxyControlBlocks(:,:,index(kk,1));
        end
        bootProxy = bootProxy(1:T_est,:);
        bootProxyControl = bootProxyControl(1:T_est,:);
        
        %center the bootstrapped residuals and proxies
        bootU = bootU - VARcentering;
        for kk = 1:size(proxy,2)
            bootProxy((bootProxy(:,kk)~=0),kk) =...
                bootProxy((bootProxy(:,kk)~=0),kk) - Proxycentering((bootProxy(:,kk)~=0),kk);
            bootProxyControl((bootProxyControl(:,kk)~=0),kk) =...
                bootProxyControl((bootProxyControl(:,kk)~=0),kk) - ProxyControlcentering((bootProxyControl(:,kk)~=0),kk);
        end
        
        % adjust for identification sample
        bootProxy = bootProxy(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:); 
        bootProxyControl = bootProxyControl(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:); 
        
        % treatment months
        bootindsR1 = zeros(nBlock*BlockSize,size(indsR1,2)); 
        bootindsR2 = zeros(nBlock*BlockSize,size(indsR2,2)); 
        for kk = 1:nBlock
            bootindsR1(1+BlockSize*(kk-1):BlockSize*kk,:) = indsR1Blocks(:,:,index(kk,1));
            bootindsR2(1+BlockSize*(kk-1):BlockSize*kk,:) = indsR2Blocks(:,:,index(kk,1));
        end
        bootindsR1 = bootindsR1(1:T_est,:);
        bootindsR2 = bootindsR2(1:T_est,:);
        
        % adjust for identification sample
        bootindsR1 = logical(bootindsR1(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:)); 
        bootindsR2 = logical(bootindsR2(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:)); 
        
        % count the number proxy variables not censored to zero
        ProxyCount(j,:) = sum(abs(bootProxy) > 0,1);
        
        if ProxyCount(j,:)<15
            continue
        end
        
        % simulate VAR starting from initial values
        Xexo = varEst.Xexo;
        
        bootData = zeros(T_est+p,nvar); 
        bootData(1:p,:) = data(1:p,:); %initial values of y, same for all j
        for i = p+1:T_est+p
            bootData(i,:)= varEst.B*[Xexo(i-p,:)'; vec(fliplr(bootData(i-p:i-1,:)'))] ...
                             + bootU(i-p,:)'; % bootstrap
        end
        
        % re-estimate the VAR
        bootvarEst = varxest(bootData,dataExo,p);
        
        % only use proxy sample for identification
        bootU = bootvarEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation  
        bootSigma = bootU'*bootU/(T-p*nvar-nexo);
    end
        
    % structural impact matrix
    bootT_OPEC = size(bootProxy(bootindsR1,:),1);
    bootT_Control = size(bootProxyControl(bootindsR2,:),1);

    bootXrIV = [(bootProxy(bootindsR1,:)-mean(bootProxy(bootindsR1,:)))/sqrt(bootT_OPEC); (bootProxyControl(bootindsR2,:)-mean(bootProxyControl(bootindsR2,:)))/sqrt(bootT_Control)];
    bootZrIV = [(bootProxy(bootindsR1,:)-mean(bootProxy(bootindsR1,:)))/sqrt(bootT_OPEC); -(bootProxyControl(bootindsR2,:)-mean(bootProxyControl(bootindsR2,:)))/sqrt(bootT_Control)];
    bootyiIV = [(bootU(bootindsR1,:)-mean(bootU(bootindsR1,:)))/sqrt(bootT_OPEC); (bootU(bootindsR2,:)-mean(bootU(bootindsR2,:)))/sqrt(bootT_Control)];

    % first stage
    bootOlsEst = olsest(bootZrIV,bootXrIV,true,true);
    bootuhat = bootOlsEst.yhat;

    % second stage
    bootb21ib11_2SLS = [bootuhat]\bootyiIV;  
    bootb1 = bootb21ib11_2SLS';      % 2 SLS coefficients

    % compute IRFs
    bootIRFs(:,:,j)  = varirfsingle(bootvarEst.B(:,1+nexo:end),bootb1,p,horizon,diffInd);
    if strcmp(shockType,'custom')
        bootIRFs(:,:,j) = bootIRFs(:,:,j)./bootIRFs(1,1,j)*shockSize;
    end 
    bootb1s(:,j) = bootb1;
    bootBs(:,:,j) = bootvarEst.B;
    bootSigma1 = bootU(bootindsR1,:)'*bootU(bootindsR1,:)/(sum(bootindsR1)); % -p*nvar-nexo
    bootShocks(:,j) = (bootb1'*inv2(bootSigma1)*bootU')'*inv2(bootb1'*inv2(bootSigma1)*bootb1);
    bootDatas(:,:,j) = bootData;
  
    j = j+1;
end

IRFsmed   = quantile(bootIRFs, 0.5, 3);
IRFslower_pe = quantile(bootIRFs, 1-alpha/2, 3)-IRFsmed+IRFs_pe;  % correct for small-sample bias
IRFsupper_pe = quantile(bootIRFs, alpha/2, 3)-IRFsmed+IRFs_pe;
IRFslower2_pe = quantile(bootIRFs, 1-alpha2/2, 3)-IRFsmed+IRFs_pe;  % correct for small-sample bias
IRFsupper2_pe = quantile(bootIRFs, alpha2/2, 3)-IRFsmed+IRFs_pe;


% get shock
Sigma1 = U(indsR1,:)'*U(indsR1,:)/(sum(indsR1)); % -p*nvar-nexo
Sigma2 = U(indsR2,:)'*U(indsR2,:)/(sum(indsR2)); % -p*nvar-nexo
Sigma = U'*U/(T-p*nvar-nexo);

oilSupplyNewsShock1 = (b1'*inv2(Sigma1)*U')'*inv2(b1'*inv2(Sigma1)*b1); 
oilSupplyNewsShock2 = (b1'*inv2(Sigma2)*U')'*inv2(b1'*inv2(Sigma2)*b1); 

oilSupplyNewsShock = oilSupplyNewsShock1; % can also use oilSupplyNewsShock2
