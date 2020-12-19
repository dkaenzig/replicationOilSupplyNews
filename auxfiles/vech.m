function [Y,H,K]=vech(X)
%
% Luca Benati
% Bank of England
% Monetary Assessment and Strategy Division
% May 2006
%
% function [Y,H,K]=vech(X)
% This program implements the vech operator for a square matrix X.
%
%                                   Input:
% X = a square matrix
%                                   Output:
% Y   = vech(X)
% H,K = the two indices corresponding to the row and column of X(H,K) for
%       each entry of Y
[R,C]=size(X);
if R~=C
    disp('Matrix is not square: program is terminated')
    return
else
    Y=[];
    H=[];
    K=[];
    for kk=1:C
        for hh=kk:R
            H=[H; hh];
            K=[K; kk];
        end
        Y=[Y; vec(X(kk:R,kk))];
    end
end