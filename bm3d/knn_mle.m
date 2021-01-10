function [ u ] = knn_mle(noisy, srch_r, sim_r, k, dist, usedist)
%knn_ml: k-nearst neighborhood maximum likelihood estimator
%
%   Poisson nonlocal maximum likelihood estimator. Use the center pixels of 
%   the k-most similar patches as samples for the maximum likelihood 
%   estimator.
%
%parameters:
%   noisy : noisy image
%   srch_r : search range
%   sim_r : similarity range
%   k : sample size
%   dist: distance used
% 
%returns:
%   u : maximum likelihood estimates 
%
%   Andr� Bindilatti
%   andre.bindilatti@dc.ufscar.br / andre_a_a@yahoo.com.br
%   

% image dimensions 
[m, n] = size(noisy);
% distance kernel
% g = ones(2*sim_r + 1);
% g = disk_kernel(sim_r);
g = fspecial('disk', sim_r);
% pad noisy image
v = padarray(noisy,[srch_r srch_r],'symmetric');
% search area
srch_area = (2*srch_r + 1)^2;
% matrix of distances
M = zeros(m, n, srch_area);
% matrix of samples
Z = zeros(m, n, srch_area);
% counter 
i = 0;

for dx = -srch_r:srch_r
for dy = -srch_r:srch_r
    
    % increment counter
    i = i + 1;
    % reference image 
    v1 = v(srch_r + (1:m), srch_r + (1:n));
    % shifted image
    v2 = v(srch_r + (1:m) + dx, srch_r + (1:n) + dy);
    % distance
    if (usedist) 
        if (dist == "euclidian_distance")
            dv = (v1 - v2).^2; %Euclidean
        end
        if (dist == "kullback_leibler_distance")
            dv = kldist(v1, v2);
        end
        if (dist == "hellinger_distance")
            dv = hellingerdist(v1, v2);
        end
        if (dist == "bhattacharryya_distance")
            dv = bhattacharyyadist(v1, v2);
        end
        if (dist == "renyi_distance")
            dv = renyidist(v1, v2, 0);
        end
    else
        dv = (v1 - v2).^2;
    end
    %dv = ((v1 - v2) * log(v1 / v2)) / 2; %Kullback
    %dv = (1-exp(-1/2*(v1+v2)+sqrt(v1.*v2))); % Hellinger
    %dv = (1/2*(v1+v2)-sqrt(v1.*v2)); %Bat
    
    d = imfilter(dv,g,'symmetric');
    % set distance matrix
    M(:,:,i) = d;
    % set samples matrix
    Z(:,:,i) = v2;
    
end
end

% output image
u = zeros(m,n);

for i = 1:m
for j = 1:n
    
    % sorted list of distances
    [~, L] = sort( M(i,j,:) );
    % nonlocal neighborhood of (i,j)
    S = Z( i, j, L(1:k) );
    % maximum likelihood estimator
    u(i,j) = sum(S)/k;
    
end
end

end

