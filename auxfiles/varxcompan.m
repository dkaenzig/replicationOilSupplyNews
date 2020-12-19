function output = varxcompan(varEst,spars)
%VARCOMPAN generates the matrices of the companion form of the VAR(p)
% build matrix F (sparse)
%
% Inputs: -varEst: estimation results from VARX(p). Relevant components:
%                  -B: n x np matrix of VAR coeffs (n x np+1 if constant is included)
%                  -Sigma: n x n covariance matrix of residuals
%                  -const: Boolean whether constant is included
%         -sparse: Boolean whether to use sparse matrices to save memory
% Output: -C: np x 1 vector of constants of companion form
%         -F: np x np matrix of VAR coeffs (without constant) of companion form
%         -VAR: np x np covariance matrix of residuals of companion form
%
% Diego R. Känzig. This version: 22/02/2018

if nargin < 2
    spars = false; % default: not sparse
end

% get inputs from varEst
B = varEst.B;
Sigma = varEst.Sigma;
nexo = varEst.nexo;

c = B(:,1:nexo);
B = B(:,nexo+1:end);

[n, pn]=size(B);

if spars
    C = spalloc(pn,nexo,pn-n*nexo);
    C(1:n,1:nexo) = c;

    F=spalloc(pn,pn,(n+1)*pn);
    F(1:n,:)=B;
    F(n+1:pn,1:pn-n)=speye(pn-n);
    
    VAR = spalloc(pn,pn,pn^2-n^2);
    VAR(1:n,1:n) = Sigma;
else
    C = zeros(pn,nexo);
    C(1:n,1:nexo) = c;

    F = zeros(pn,pn);
    F(1:n,:) = B;
    F(n+1:pn,1:pn-n) = eye(pn-n);
    
    VAR = zeros(pn,pn);
    VAR(1:n,1:n) = Sigma;
end

output.C = C;
output.F = F;
output.VAR = VAR;

end
