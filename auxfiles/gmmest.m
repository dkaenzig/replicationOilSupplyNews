function output = gmmest(X,Z,y,const,robust,lag)
%% This function implements efficient GMM of y on variables in X using variables Z as instruments
% Inputs: -X: matrix containing the independent variables (stacked as column
%             vectors)
%         -Z: matrix conaining the instruments (stacked as column vectors)
%         -y: column vector of the dependent variable
%         -const: a boolean indicating whether to include a constant or not
%         -robust: indicator for:
%                  -standard SE (assuming conditional homoskedasticity) (0),
%                  -White (heteroskedasticity robust) SE (1) or 
%                  -HAC (heteroskedasticity and autocorrelation consistent) SE (2)
%                   Note: HAC compted using Newey-West estimator and
%                   bartlett window
%         -lag: number of lags used in autocorrelation structure (only
%               needed if robust = 2)
% Outputs:-output: struct containing the fields
%                   -bhat: vector containing OLS estimates
%                   -yhat: fitted values
%                   -uhat: residuals
%                   -rmse: root mean squared error
%                   -n: sample size
%                   -k: number of regressors   
%                   -l: number of instruments 
%                   -varbhat: variance-covariance matrix of bhat (depending
%                             on the option for robust chosen)
%                   -lag: number of lags used in autocorrelation structure
%
% Diego R. Känzig. This version: 21/02/2018

    n = size(X,1);
    if nargin==5 && robust==2
        lag = round(4*(n/100)^(2/9));  % set it to Newey and West rule of thumb
    elseif nargin==5 % && robust~=2
        lag = 0;
    end
    if nargin==4
        robust = 0;
        lag = 0;
    end
    if nargin==3
        robust = 0;
        lag = 0;
        const = 0;
    end

    if const
        X = [X ones(n,1)];
        Z = [Z ones(n,1)];
    end
    k = size(X,2);
    l = size(Z,2);
    
    % 1: get residuals from 2SLS
    PZ = Z/(Z'*Z)*Z';
    b2sls = (X'*PZ*X)\X'*PZ*y;   % consistent estimator
    uhat2sls = y - X*b2sls;
    
    % 2: compute efficient weighting matrix
    if robust==0
        Shat = 1/(n-k)*(uhat2sls'*uhat2sls)*(Z'*Z); 
    elseif robust==1
        Shat = zeros(k);      % more efficient with loop
        for ii = 1:n
            Shat = Shat+uhat2sls(ii)^2*Z(ii,:)'*Z(ii,:);
        end
        Shat = n/(n-k)*Shat;
    elseif robust==2
        Shat0 = zeros(k);      % more efficient with loop   
        for ii = 1:n
            Shat0 = Shat0+uhat2sls(ii)^2*Z(ii,:)'*Z(ii,:);
        end
        Shat = Shat0;
        for il = 1:lag
            Shatl = zeros(k); 
            for ii = 1+il:n
                Shatl = Shatl+(1-il/(lag+1))*uhat2sls(ii)*uhat2sls(ii-il)*(Z(ii,:)'*Z(ii-il,:)+Z(ii-il,:)'*Z(ii,:)); % Bartlett window
            end
            Shat = Shat + Shatl;
        end
        Shat = n/(n-k)*Shat;
    end
    
    % efficient weighting matrix
    W = inv2(Shat); %inv(Shat);
    
    % gmm estimator
    bhat = (X'*Z*W*Z'*X)\(X'*Z*W*Z'*y);
    yhat  = X*bhat;
    uhat  = y-yhat;
    rmse  = sqrt(1/(n-k)*(uhat'*uhat));
    
    varbhat = inv2(X'*Z*W*Z'*X);
    
    output.bhat    = bhat;
    output.yhat    = yhat;
    output.uhat    = uhat;
    output.rmse    = rmse;
    output.n       = n;
    output.k       = k;
    output.l       = l;
    output.varbhat = varbhat;
    output.lag     = lag;
     
end
