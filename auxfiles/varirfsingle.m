function [IRF] = varirfsingle(B,invA0i,p,hor,cum)
% computes IRFs of a VAR model of the form
% Y_t=C+B*Y_{t-1}+A0*V_t
% for a prespecified horizon
%
% Inputs: -B: Coefficient matrix of the VAR of dimension n x (n*p) (exclude deterministics)
%         -invA0i: i th column of the structural impact matrix
%         -hor: Horizon
%
% Output: -IRF: The IRFs for the variables in Y_t with respect to the i th
%          shock stacked in columns
%
% Diego R. Känzig. This version: 10/02/2018

if nargin < 5
    cum = zeros(size(invA0i,1),1);
end
    
IRF = zeros(size(invA0i,1),p+hor+1);
IRF(:,p+1) = invA0i;
for tt = p+2:p+hor+1
    IRF(:,tt) = B*vec(fliplr(IRF(:,tt-p:tt-1)));
end
%
IRF=IRF(:,p+1:end)';

for i = find(cum) % compute CIRF if requested
    IRF(:,i) = cumsum(IRF(:,i));
end
