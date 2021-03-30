function rmse=calrmse(x,y)
error=(x-y).^2/size(x,1)/size(x,2);
rmse=sqrt(sum(error(:)));