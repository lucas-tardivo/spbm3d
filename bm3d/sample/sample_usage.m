% Default parameters
peak = 100;
theta = 7;
lambdaHard = 0.5;
rng(561996);

% Lambda parameters
srch_r = 7;
sim_r = 5;
k = 30;

% Read a mat image
% IMG = load('Carcinoma_Angular_BUS_1.mat');
% I = IMG.Data.Full.EchoModel;

% Read a png image
 I = imread('baboon.png'); 

% Prepare image
noiseless = double(I);
[x,y] = size(noiseless);
xstart = round((x-256)/2);
ystart = round((y-256)/2);
noiseless = noiseless(xstart:(xstart+(256-1)), ystart:(ystart+(256-1)));

% Simulate noise
[Q,noisyimage] = poissnoise(im2double(noiseless), peak); 
noisy = noisyimage;

% Set Step 1 Parameters
step1Parms.use_std = 1;
step1Parms.transform = 1;
step1Parms.sd_thr = 0.1;
step1Parms.s = 1;
step1Parms.b = 1.1;
step1Parms.lambdaHard3D = lambdaHard;
step1Parms.tauMatch = 0;

step1Parms.distance_func = 'hellinger_distance'; 
% -- Distance functions --
% euclidean_distance 
% bhattacharryya_distance 
% kullback_leibler_distance
% hellinger_distance
% renyi_distance

step1Parms.mle_func = 'he_greenshields_sigma_mle';
step1Parms.blk_func = 'sd';
step2Parms = step1Parms; % Same parameters for step 2

% Estimate lambda
u = knn_mle(noisy, srch_r, sim_r, k, step1Parms.distance_func, false);

% Stochastic BM3D
filtered_sd = sd_bm3d(noisy, theta*std(std(noisy)), step1Parms, step2Parms, u);
filtered_psnr = psnr(noiseless, filtered_sd, 255, Q);
filtered_ssim = ssim(filtered_sd, noiseless);

% Show results
fprintf(strcat(" \n Stochastic Result: ", num2str(filtered_psnr)), "db");
fprintf(strcat(" \n Stochastic SSIM: ", num2str(filtered_ssim)));

figure;
subplot(2,2,1);imagesc(noiseless);colormap(gray(255));title(strcat('Noiseless(ssim=', num2str(ssim(noiseless, noiseless)), ')'));
subplot(2,2,2);imagesc(noisy);colormap(gray(255));title(strcat('Noisy(ssim=', num2str(ssim(noisy, noiseless)), ')'));
subplot(2,2,4);imagesc(filtered_sd);colormap(gray(255));title(strcat('Stochastic PSNR=', num2str(filtered_psnr),'dB)'));

