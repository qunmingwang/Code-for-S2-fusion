function [RMSE,CC]=evaluate_relation(realdata,predictdata);
[a,b]=size(realdata);
RMSE=sum(sum((realdata-predictdata).^2));RMSE=sqrt(RMSE/(a*b));
C_1=sum(sum(predictdata.*realdata))-a*b*mean(mean(predictdata))*mean(mean(realdata));
C_2=sum(sum(predictdata.^2))-a*b*mean(mean(predictdata))^2;
C_3=sum(sum(realdata.^2))-a*b*mean(mean(realdata))^2;
CC=C_1/sqrt(C_2*C_3);