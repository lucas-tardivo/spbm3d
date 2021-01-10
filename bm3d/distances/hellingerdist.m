function [ d ] = hellingerdist( X, Y )

d = (1-exp(-1/2*(X+Y)+sqrt(X.*Y)));

end

