function LogitX=Logit(X)
if any(vec(squeeze(any(X)<0)))==1
    disp('Error: at least one entry is smaller than zero')
    LogitX=-9999;
    return
end
if any(vec(squeeze(any(X)>1)))==1
    disp('Error: at least one entry is greater than 1')
    LogitX=-9999;
    return
end
%
LogitX=log(X)-log(1-X);
%



