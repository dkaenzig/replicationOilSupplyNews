function qdate = mtoqdate(mdate)
%MTOQDATE Transforms monthly date ('YYYYMXX') to quarterly date ('YYYYQX')

year = mdate(1:4);
if strcmp(mdate(end-1:end),'01') || strcmp(mdate(end-1:end),'02') || strcmp(mdate(end-1:end),'03')
    quarter = '1';
elseif strcmp(mdate(end-1:end),'04') || strcmp(mdate(end-1:end),'05') || strcmp(mdate(end-1:end),'06')
    quarter = '2';
elseif strcmp(mdate(end-1:end),'07') || strcmp(mdate(end-1:end),'08') || strcmp(mdate(end-1:end),'09')
    quarter = '3';
elseif strcmp(mdate(end-1:end),'10') || strcmp(mdate(end-1:end),'11') || strcmp(mdate(end-1:end),'12')
    quarter = '4';
end

qdate = strcat(year,'Q',quarter);

end

