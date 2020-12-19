function output = ardlest(y,x,p,q,const)
%ARDLEST estimates an ARDL(p,q) model for time series y and exogenous series
% in x (different time series stacked in columns using OLS
%
% Inputs: -y: engogenous time series
%         -x: exogenous time series (stacked in columns) [x1 x2 ...]
%         -p: lag order of endog time series
%         -q: lag order of exogenous time series
%         -const: Boolean whether to include a constant
% Outputs:-output: struct that contains all relevant output:
%                  -nobs: number of obs used in estimation
%                  -dfm: degrees of freedom of the model
%                  -bhat: the estimated coefficients [b_endo; b_exo1;
%                                                     b_exo2; ...; b_exon; const]
%                  -yhat: the fitted values
%                  -resid: the residuals
%                  -ssr: sum of squared residuals
%                  -rmse: root mean squared error
%                  -varbhat: the covariance matrix of the coefficients
%                  -R2: the R squared
%                  -adjR2: the adjusted R squared
%                  -logL: the log likelihood of the model
%                  -AIC, BIC, HQ: different information criteria for lag
%                                 length selection
%
% Diego R. Känzig. This version: 26/04/2017

if nargin<5
    const = true;
end

order = max(p,q);
X = nan(size(y,1)-order,p+q);
for i = 1:p
    X(:,i)=y(1+order-i:end-i);
end
for j = 1:size(x,2)
    for i = 1:q
        X(:,p*j+i)=x(1+order-i:end-i,j);
    end
end
T = size(X,1);
if const 
    X = [X ones(T,1)];
end
k = size(X,2);
y = y(1+order:end);

invXX = inv(X'*X);
bhat = invXX*(X'*y);
yhat = X*bhat;
resid = y-yhat;
ssr = resid'*resid;
rmse = sqrt(1/(T-k)*ssr);  
varbhat = rmse^2*invXX; 
R2 = 1-ssr/((y-mean(y))'*(y-mean(y))); 
adjR2 = 1-(ssr/(T-k))/((y-mean(y))'*(y-mean(y))/(T-1));
logL = -T/2*log(2*pi/T)-T/2-T/2*log(ssr);
AIC = -2*logL+k*2;
BIC = -2*logL+k*log(T);
HQ = -2*logL+k*2*log(log(T));

% assign outputs
output.nobs = T;
output.dfm = p+q;
output.bhat = bhat;
output.yhat = yhat;
output.resid = resid;
output.ssr = ssr;
output.rmse = rmse; 
output.varbhat = varbhat;
output.R2 = R2;
output.adjR2 = adjR2;
output.logL = logL;
output.AIC = AIC;
output.BIC = BIC;
output.HQ = HQ;

end

