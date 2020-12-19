%% Run the Proxy VAR

% robustness: check for optimal lag length according to information criteria
if doLagCrit
    if strcmp(dataFrequency,'M')
        pmax = 36;
    elseif strcmp(dataFrequency,'Q')
        pmax = 12;
    end
    infoCrit = nan(pmax,3);
    for ii = 1:pmax
        varSel = varxest(data,dataExo,ii,false);      % don't do dfm correction for Sigma
        infoCrit(ii,:) = [varSel.AIC varSel.BIC varSel.HQ];
    end
    [~,pSel] = min(infoCrit);     % these lag orders can be used as a robustness check

    if strcmp(p,'aic')
        p = pSel(1);
    end
end

% run reduced-form VAR 
varEst = varxest(data,dataExo,p);

% identification using the covariance structure between proxy and
% reduced-form residuals a la Mertens and Ravn (2013)

% only use proxy sample for identification (potentially a subset of the
% estimation sample)
nexo = size(dataExo,2);
U = varEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation  
Sigma = U'*U/(T-p*nvar-nexo);

% run first stage (k variables to be instrumented with size(proxy,2)
% proxies. Need that size(proxy,2)>=k)
uhat = [];
for ik = 1:k
    inam = strcat('r',num2str(ik));
    olsEst.(inam) = olsest(proxy,U(:,ik),true,true);

    disp('First-stage:')
    fprintf('F-stat: %4.3f, p-value: %4.3f, F-stat (robust): %4.3f, p-value: %4.3f, R^2: %4.3f, R^2 (adj): %4.3f \n',olsEst.(inam).F,olsEst.(inam).Fpval,olsEst.(inam).Frobust,olsEst.(inam).Frobustpval,olsEst.(inam).R2,olsEst.(inam).R2adj)
    uhat = [uhat olsEst.(inam).yhat];
end

% second stage
b21ib11_2SLS    =   [ones(length(proxy),1) uhat]\U(:,k+1:end);  
b21ib11 = b21ib11_2SLS(2:end,:)';      % 2 SLS coefficients               
Sig11   = Sigma(1:k,1:k);
Sig21   = Sigma(k+1:nvar,1:k);
Sig22   = Sigma(k+1:nvar,k+1:nvar);
ZZp     = b21ib11*Sig11*b21ib11'-(Sig21*b21ib11'+b21ib11*Sig21')+Sig22;
b12b12p = (Sig21- b21ib11*Sig11)'*(ZZp\(Sig21- b21ib11*Sig11));
b11b11p = Sig11-b12b12p;
b11     = sqrt(b11b11p);
b1      = [b11; b21ib11*b11];
b1unit  = [1; b21ib11]*shockSize;

% Reliability (alternative to F stat)
Sigmm   = proxy'*proxy/T;
ED      = eye(k)*sum(sum(proxy,2)~=0)/T;
mu1     = proxy'*U(:,1:k)/T;
PhiPhip = mu1*inv2(b11b11p)*mu1';
RM      = inv2(Sigmm)*PhiPhip*inv2(ED);
RMeigs  = sort(eig(RM));

% Alternative measure if k = 1 (this measure is preferred by Mertens and Ravn)
if k==1 && np==1
    Bi = [1/b11+(Sig21-Sig11*b21ib11)'*inv2(ZZp)/b11*b21ib11 -(Sig21-Sig11*b21ib11)'*inv2(ZZp)/b11];
    et = (Bi*U')';

    PHI = mu1/b11;
    GAM = inv2(ED)*PHI;
    E  = GAM*et(sum(proxy,2)~=0);
    V  = proxy(sum(proxy,2)~=0)-E;
    RM_2 = inv2(E'*E+V'*V)*E'*E;
end

% Lunsford alternative to Gertler Karadi F-stat
olsEstLunds = olsest(U,proxy,true,true);

% compute IRFs to shock
IRFs_pe = varirfsingle(varEst.B(:,1+nexo:end),b1,p,horizon,diffInd);
if strcmp(shockType,'custom')
    IRFs_pe = IRFs_pe./IRFs_pe(1,1)*shockSize;
end

% compute FEVDs to shock
if getFEVDs
    varEst.Sigma = Sigma;   % use Sigma based on truncated sample for FEVD ?
    FEVDs_pe = varxfevdsingle(varEst,b1,horizon+1); 
end

% compute the confidence bands using bootstrapping
bootIRFs = nan(horizon+1,nvar,nsim);
bootFEVDs = nan(horizon+1,nvar,nsim);
bootFstat = nan(nsim,2);
bootBs = nan(nvar,size(varEst.B,2),nsim);
bootb1s = nan(nvar,nsim);
bootShocks = nan(T,nsim);
bootDatas = zeros(T+p,nvar,nsim); 
ProxyCount = zeros(nsim,size(proxy,2));
T_est = varEst.T; % length of estimation sample

if strcmp(bootType,'mbb1block')
    % if identification sample is shorter that estimation sample, censor
    % unobserved values to zero
    proxyLong = zeros(T_est, size(proxy,2));
    proxyLong(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:) =  proxy;

    BlockSize = round(5.03*T_est^0.25);
    nBlock = ceil(T_est/BlockSize);
    VARBlocks = zeros(BlockSize,nvar,T_est-BlockSize+1);
    ProxyBlocks = zeros(BlockSize,size(proxyLong,2),T_est-BlockSize+1);
    for j = 1:T_est-BlockSize+1
        VARBlocks(:,:,j) = varEst.U(j:BlockSize+j-1,:);
        ProxyBlocks(:,:,j) = proxyLong(j:BlockSize+j-1,:);
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
    for j = 1:BlockSize 
        subProxy = proxyLong(j:T_est-BlockSize+j,:);
        %Proxycentering(j,:) = mean(subProxy((subProxy(:,1) ~= 0),1),1); 
        % account for non-zero mean instrument:
        Proxycentering(j,:) = mean(subProxy((subProxy(:,1) ~= 0),1),1) - mean(proxyLong((proxyLong(:,1) ~= 0),1),1);
    end
    Proxycentering = repmat(Proxycentering,[nBlock,1]);
    Proxycentering = Proxycentering(1:T_est,:);
    
end

j = 1;
while j <= nsim
    % generate artificial data
    
    if strcmp(bootType,'wild')
        % Wild bootstrap using Rademacher (should not use because invalid)
        
        rr = 1-2*(rand(T_est,1)>0.5);       % use wild bootstrap based on Rademacher distribution
        bootU = (varEst.U).*(rr*ones(1,nvar)); % draw from reduced-form shocks
        bootProxy = proxy.*(rr(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:)*ones(1,size(proxy,2))); % draw from proxy (note the adjustment because of shorter sample)

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
        % and the proxy
        
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
        bootU = bootUhat(index2,:);
        bootSigma = bootU'*bootU/(T-p*nvar-nexo);
        
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
        for kk = 1:nBlock
            bootProxy(1+BlockSize*(kk-1):BlockSize*kk,:) = ProxyBlocks(:,:,index(kk,1));
        end
        bootProxy = bootProxy(1:T_est,:);

        %center the bootstrapped residuals and proxies
        bootU = bootU - VARcentering;
        for kk = 1:size(proxy,2)
            bootProxy((bootProxy(:,kk)~=0),kk) =...
                bootProxy((bootProxy(:,kk)~=0),kk) - Proxycentering((bootProxy(:,kk)~=0),kk);
        end
        
        % adjust for identification sample
        bootProxy = bootProxy(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:); 
        
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
    % 2SLS of U1 on U2 using proxy as an instrument (include a constant for the case that proxy is not demeaned)

    % first stage
    bootOlsEst = olsest(bootProxy,bootU(:,1),true,true);
    bootuhat = bootOlsEst.yhat;

    bootFstat(j,:) = [bootOlsEst.F bootOlsEst.Frobust];
    
    % second stage
    bootb21ib11_2SLS    =   [ones(length(proxy),1) bootuhat]\bootU(:,k+1:end);  
    bootb21ib11 = bootb21ib11_2SLS(2:end,:)';      % 2 SLS coefficients    

    bootSig11   = bootSigma(1:k,1:k);
    bootSig21   = bootSigma(k+1:nvar,1:k);
    bootSig22   = bootSigma(k+1:nvar,k+1:nvar);
    bootZZp     = bootb21ib11*bootSig11*bootb21ib11'-(bootSig21*bootb21ib11'+bootb21ib11*bootSig21')+bootSig22;
    bootb12b12p = (bootSig21- bootb21ib11*bootSig11)'*(bootZZp\(bootSig21- bootb21ib11*bootSig11));
    bootb11b11p = bootSig11-bootb12b12p;
    bootb11     = sqrt(bootb11b11p);
    bootb1      = [bootb11; bootb21ib11*bootb11];

    % compute IRFs
    bootIRFs(:,:,j)  = varirfsingle(bootvarEst.B(:,1+nexo:end),bootb1,p,horizon,diffInd);
    if strcmp(shockType,'custom')
        bootIRFs(:,:,j) = bootIRFs(:,:,j)./bootIRFs(1,1,j)*shockSize;
    end 
    if getFEVDs
        bootvarEst.Sigma = bootSigma;
        bootFEVDs(:,:,j) = varxfevdsingle(bootvarEst,bootb1,horizon+1); 
    end
    bootBs(:,:,j) = bootvarEst.B;
    bootb1s(:,j) = bootb1;
    bootShocks(:,j) = (bootb1'*inv2(bootSigma)*bootU')';
    bootDatas(:,:,j) = bootData;
    
    j = j+1;
end

if strcmp(bootType,'mbb1block') || strcmp(bootType,'iid')
    % display the fewest non-censored proxy observations in the bootstrap
    disp('Minimum non-censored proxy observations in bootstrap:')
    disp(min(ProxyCount))
end

% rescale bootstrapped IRFs (center around sample estimates
% or perform bias correction) and get quantiles
IRFsmed_pe = quantile(bootIRFs, 0.5, 3);

IRFsupper_pe = quantile(bootIRFs, 1-alpha/2, 3)-IRFsmed_pe+IRFs_pe;  % rescaling does not affect the ordering, thus it is fine to take quantile first
IRFslower_pe = quantile(bootIRFs, alpha/2, 3)-IRFsmed_pe+IRFs_pe;

IRFsupper2_pe = quantile(bootIRFs, 1-alpha2/2, 3)-IRFsmed_pe+IRFs_pe;  % rescaling does not affect the ordering, thus it is fine to take quantile first
IRFslower2_pe = quantile(bootIRFs, alpha2/2, 3)-IRFsmed_pe+IRFs_pe;

% rescale bootstrapped FEVDs (a bit more subtle because we have to
% make sure that corrected FEVD lie within [0,1]. Use logistic
% functions) and get quantiles

if getFEVDs
    FEVDslogit = Logit(FEVDs_pe);
    bootFEVDslogit = Logit(bootFEVDs);
    FEVDsmedlogit = quantile(bootFEVDslogit, 0.5, 3);

    FEVDsupper_pe = InverseLogit(quantile(bootFEVDslogit, 1-alpha/2, 3)-FEVDsmedlogit+FEVDslogit); 
    FEVDslower_pe = InverseLogit(quantile(bootFEVDslogit, alpha/2, 3)-FEVDsmedlogit+FEVDslogit);  
end

% get shock as in Stock and Watson (2018)
if strcmp(shockType,'custom')
    % unit normalization
    oilSupplyNewsShock = (b1unit'*inv2(Sigma)*U')'*inv2(b1unit'*inv2(Sigma)*b1unit);
else
    % one sd normalization
    oilSupplyNewsShock = (b1'*inv2(Sigma)*U')'*1; 
end

