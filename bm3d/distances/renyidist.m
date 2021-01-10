function [ d ] = renyidist( X, Y, beta )

 b = beta;

 exp1=exp(power(Y,1-b).*power(X,b)-((1-b)*Y+b*X));
 exp2=exp(power(X,1-b).*power(Y,b)-((1-b)*X+b*Y));
 
 d = (1/(b-1)*log((exp1+exp2)/2));
 
 end

