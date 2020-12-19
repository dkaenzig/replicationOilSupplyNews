function output = olsest(X,y,const,robust,lag)
%% This function implements OLS of y on variables in X
% Inputs: -X: matrix containing the independent variables (stacked as column
%             vectors)
%         -y: column vector of the dependent variable
%         -const: a boolean indicating whether to include a constant or not
%         -robust: indicator for standard SE (assuming conditional homoskedasticity) (0),
%                  White (heteroskedasticity robust) SE (1) or 
%                  HAC (heteroskedasticity and autocorrelation consistent) SE (2)
%                  Note: HAC compted using Newey-West estimator and
%                  bartlett window
%         -lag: number of lags used in autocorrelation structure (only
%               needed if robust = 2)
% Outputs:-output: struct containing the fields
%                   -bhat: vector containing OLS estimates
%                   -yhat: fitted values
%                   -resid: residuals
%                   -SSR: sum of squared residuals
%                   -SST: sum of squares total
%                   -SSE: sum of squares explained
%                   -rmse: root mean squared error
%                   -R2: R squared
%                   -R2: adj R squared
%                   -n: sample size
%                   -k: number of regressors              
%                   -varbhat: variance-covariance matrix of bhat
%                   -varbhatrobust: robust (White or HAC) variance-covariance matrix of bhat
%                   -F: F statistic for test of overall significance (only
%                       well defined if there is a constant
%                   -Fpval: p-value of overall significance test
%                   -Frobust: robust F statistic (based on large sample OLS distribution)
%                   -Frobustpval: p-value of robust F statistic
%                   -lag: number of lags used in autocorrelation structure
%
% Diego R. Känzig. This version: 21/02/2018

    n = size(X,1);
    if nargin==4 
        if robust==2
            lag = round(4*(n/100)^(2/9));  % set it to Newey and West rule of thumb
        elseif robust~=2  
            lag = 0;
        end
    end
    if nargin==3
        robust = 0;
        lag = 0;
    end
    if nargin==2
        robust = 0;
        lag = 0;
        const = 0;
    end

    if const
        X = [X ones(n,1)];
    end
    k = size(X,2);
    invXpX = (X'*X)\eye(k); % inv(X'*X)
    bhat  = invXpX*(X'*y);
    yhat  = X*bhat;
    resid = y-yhat;
    SSR   = resid'*resid;
    SST   = (y-mean(y))'*(y-mean(y));
    SSE   = (yhat-mean(y))'*(yhat-mean(y));
    rmse  = sqrt(1/(n-k)*SSR);
    R2    = 1 - SSR/SST; 
    R2adj = 1 - (n-1)/(n-k)*SSR/SST;
    
    varbhat = rmse^2*invXpX; 
    if robust==1
        Shat = zeros(k);      % more efficient with loop
        for ii = 1:n
            Shat = Shat+resid(ii)^2*X(ii,:)'*X(ii,:);
        end
        varbhatrobust = n/(n-k)*invXpX*Shat*invXpX;
        %varbhat= n/(n-k)*invXX*(X'*diag(resid.^2)*X)*invXX;
    elseif robust==2
        Shat0 = zeros(k);      % more efficient with loop   
        for ii = 1:n
            Shat0 = Shat0+resid(ii)^2*X(ii,:)'*X(ii,:);
        end
        Shat = Shat0;
        for l = 1:lag
            Shatl = zeros(k); 
            for ii = 1+l:n
                Shatl = Shatl+(1-l/(lag+1))*resid(ii)*resid(ii-l)*(X(ii,:)'*X(ii-l,:)+X(ii-l,:)'*X(ii,:)); % Bartlett window
            end
            Shat = Shat + Shatl;
        end
        varbhatrobust = n/(n-k)*invXpX*Shat*invXpX;
    end

    % F test of overall significance
    R = zeros(k-1,k);
    R(1:k-1,1:k-1) = eye(k-1);  % all coeffs zero except constant
    q = zeros(k-1,1);

    F = 1/(k-1)*(R*bhat-q)'*inv(rmse^2*R*invXpX*R')*(R*bhat-q);         % small sample
    Fpval       = 1 - fcdf(F,k-1,n-k);
    
    % see https://stats.stackexchange.com/questions/93787/f-test-formula-under-robust-standard-error
    if robust>0
        Frobust = 1/(k-1)*(R*bhat-q)'*inv(n/(n-k)*R*invXpX*Shat*invXpX*R')*(R*bhat-q);    % large sample (Wald test)
        % Frobust = 1/(k-1)*(R*bhat-q)'*inv(R*varbhatrobust*R')*(R*bhat-q);    % large sample (Wald test)
        Frobustpval = 1 - fcdf(Frobust,k-1,n-k);
    end
    
    output.bhat    = bhat;
    output.yhat    = yhat;
    output.resid   = resid;
    output.SSR     = SSR;
    output.SST     = SST;
    output.SSE     = SSE;
    output.rmse    = rmse;
    output.R2      = R2;
    output.R2adj   = R2adj;
    output.n       = n;
    output.k       = k;
    output.varbhat = varbhat;
    if robust>0
        output.varbhatrobust = varbhatrobust;
    end
    output.F       = F;
    output.Fpval   = Fpval;
    if robust>0
        output.Frobust     = Frobust;
        output.Frobustpval = Frobustpval;
    end
    output.lag     = lag;
     
end
