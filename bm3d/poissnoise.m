function [Q,noisy] = poissnoise(img, peak)

Q = max(max(img)) / peak;   
img_low = img / Q;
img_low(img_low == 0) = min(min(img_low(img_low > 0)));
noisy = poissrnd(img_low);

end

