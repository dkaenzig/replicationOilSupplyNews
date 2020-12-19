function output = varxest(data,dataexo,p,dfmcorr)
%VARXEST estimates a VARX(p) model for the series contained in data and dataexo by OLS
% Inputs: - data: the endogeneous variables stacked in columns
%         - dataexo: the exogenous variables in the VAR (including deterministics)
%         - p: the lag order of the VAR
%         - const: boolean whether to include a constant (true vs false)
% Output: - B: OLS coefficients (constant, Phi_1, Phi_2, ... , Phi_p) for
%              the n series (n x const+np)
%         - Yhat: Fitted values
%         - U: Residuals
%         - Sigma: Variance-covariance matrix residuals
%
% Diego R. Känzig. This version: 10/00/2018

if nargin < 4
    dfmcorr = true; % default: use dfm correction
end

% Define the data matrices
Y = data(p+1:end,:);        % define left hand variables
[T,n] = size(Y);

X = lagmatrix(data,1:p);    % define right hand variables
X = X(p+1:end,:); 
Xexo = dataexo(p+1:end,:); 
[~,nexo] = size(Xexo);

X=[Xexo X];

B    = (X'*X)\(X'*Y); % OLS estimates
Yhat = X*B;         % fitted values

U = Y-Yhat;                   % residuals
if dfmcorr
    Sigma = U'*U/(T-p*n-nexo);   % covariance matrix of U (using small-sample dof correction)
else
    Sigma = U'*U/T;   % covariance matrix of U (MLE estimate)
end

VARB = kron(inv(X'*X),Sigma); % covariance matrix of B'

logL = -(T*n/2)*log(2*pi)+T/2*log(det(inv(Sigma)))-(T*n/2); % Log-likelihood

AIC = log(det(Sigma))+2*(n^2*p+nexo)/T;                    % Information criteria
BIC = log(det(Sigma))+(n^2*p+nexo)/T*log(T);
HQ  = log(det(Sigma))+2*(n^2*p+nexo)/T*log(log(T));

% define outputs
output.Y     = Y;
output.X     = X;
output.Xexo  = Xexo;
output.B     = B';
output.Yhat  = Yhat;
output.U     = U;
output.Sigma = Sigma;
output.VARB  = VARB;
output.logL  = logL;
output.AIC   = AIC;
output.BIC   = BIC;
output.HQ    = HQ;
output.T     = T;
output.n     = n;
output.p     = p;
output.nexo  = nexo;

end

