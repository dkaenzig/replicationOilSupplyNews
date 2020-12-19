function FEVDi = varxfevdsingle(varEst,invA0i,hor)
% computes the FEVD of a VAR(p) model
% Y_t=B0+B1*Y_{t-1}+...+Bp*Y_{t-p}+A0*V_t, with V_t~N(0,1)
% for a prespecified horizon
%
% Inputs: -varEst: Estimation results of the VAR. Relevant components:
%                  -B: Coefficient matrix of the VAR: B=[B1 B2 ... Bp]
%                  -nexo: Number of exogenous variables
%         -invA0i: ith column of structural impact matrix satisfying
%                  invA0*invA0' = Sigma
%         -hor: Horizon
%
% Output: -FEVD: The FEVDs for the variables in Y_t with respect to the ith
%          structural shock (hor x nvar)
%
% Reference: Benati slides 10 Oct 13
% Diego R. Känzig. This version: 22/02/2018

% define inputs in varEst
B = varEst.B;
B = B(:,varEst.nexo+1:end);
Sigma = varEst.Sigma;

% get dimension and lag order
[N,pn] = size(B);
p = pn/N;

% compute companion form
F = varcompan(B);

% perform variance decomposition
% total MSE
VARk = [Sigma zeros(N,N*(p-1)); zeros(N*(p-1),N) zeros(N*(p-1),N*(p-1))];

MSE(:,:,1) = VARk;
for kk = 2:hor
   VARk = F*VARk*F';
   MSE(:,:,kk) = MSE(:,:,kk-1) + VARk;
end;

% MSEs per shock
FEVD = zeros(N,hor);
VARk_j = [invA0i*invA0i' zeros(N,N*(p-1)); zeros(N*(p-1),N) zeros(N*(p-1),N*(p-1))];
MSE_j(:,:,1) = VARk_j;
for kk = 2:hor
    VARk_j = F*VARk_j*F';
    MSE_j(:,:,kk) = MSE_j(:,:,kk-1) + VARk_j;   
end

% Compute the Forecast Error Covariance Decomposition
FECD = MSE_j./MSE;

% Select only the variance terms
for nn = 1:hor
   FEVD(:,nn) = diag(FECD(1:N,1:N,nn));  
end

FEVDi = permute(FEVD,[2 1]);

end

% Subfunctions
function F = varcompan(B)
% build matrix F (sparse)
% Output:  F - np x np matrix
% Input: phi - n x np matrix of VAR coeffs (without constant)

    [n, pn]=size(B);
    F=spalloc(pn,pn,(n+1)*pn);
    F(1:n,:)=B;
    F(n+1:pn,1:pn-n)=speye(pn-n);
end