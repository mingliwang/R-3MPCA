function PROCv(stdGau)
%
%

for t = 1: 30
    maxNumCompThreads(1);
end
addpath(genpath('MPCAv'))
addpath(genpath('PUBLIC'))




% load mri
Im = importdata('dataset/t1_icbm_normal_1mm_pn0_rf0.mat');
Im = single(Im);


% noisy images
Nm = imNoise(Im, stdGau, 'Gaussian');

% the algorithm
tic;
[Jm, ~] = ALGOv(Im, Nm, stdGau);
tim = toc;


ind = 90;
psnr_slice = psnr(uint8(Jm(:, :, ind)), uint8(Im(:, :, ind)));
fprintf('Gaussian noise level %d, psnr of denoised %dth slice = %.4f \n', stdGau, ind, psnr_slice)


rmpath(genpath('MPCAv'))
rmpath(genpath('PUBLIC'))
end


