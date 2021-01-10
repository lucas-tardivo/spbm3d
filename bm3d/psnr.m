function [psnr]=psnr(A, B, MAX, Q)

 %A = A / Q; 
 %B = B / Q;
 %d = sum((A(:)-B(:)).^2) / numel(A);
 %psnr = 10*log10((MAX*MAX)/d);
 
 B = B * Q;
 MSE = mean2((A-B).*(A-B));
 
 psnr = 10*log10(MAX^2/MSE);

end
 