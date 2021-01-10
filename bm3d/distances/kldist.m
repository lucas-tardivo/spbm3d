function [ d ] = kldist( X, Y )

d = (1/2*(X-Y).*log(X./Y));

end

