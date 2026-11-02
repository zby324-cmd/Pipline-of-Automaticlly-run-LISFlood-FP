function [b,s2,r2]=regressNoLic(y,X)
b=X\y;
mu=X*b;
s2=(y-mu)'*(y-mu)./(length(y)-length(b));
r=corrcoef(y,mu);
r2=r(1,2)^2;
end