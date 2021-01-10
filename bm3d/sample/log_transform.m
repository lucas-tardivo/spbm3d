function [ldata] = log_transform(data)
    ldata=log(data+1); 
    ldata(isnan(ldata)|isinf(ldata)) = 0;   
end