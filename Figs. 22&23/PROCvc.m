function PROCvc(stdGau)
%
%

for t = 1: 30
    maxNumCompThreads(1);
end
addpath(genpath('MPCAvc'))
addpath(genpath('PUBLIC'))




% load video
Im = importdata('dataset/news_cif.mat');
Im = single(Im);


% noisy images
Nm = imNoise(Im, stdGau, 'Gaussian');

% the algorithm
tic;
[Jm, ~] = ALGOvc(Im, Nm, stdGau);
tim = toc;


ind = 10;                                                             % original video 19th frame
psnr_frame = psnr(uint8(Jm(:, :, :, ind)), uint8(Im(:, :, :, ind)));
fprintf('Gaussian noise level %d, psnr of denoised %dth frame = %.4f \n', stdGau, 2*ind-1, psnr_frame)


rmpath(genpath('MPCAvc'))
rmpath(genpath('PUBLIC'))
end


