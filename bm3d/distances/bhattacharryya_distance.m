function [ d ] = bhattacharryya_distance( X, Y )

 d = (1/2*(X+Y)-sqrt(X.*Y));
 
end

