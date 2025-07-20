function X_hat3 = PRMc_tvs(X_orig, rng_s)
%
%

for t = 1: 30
    maxNumCompThreads(1);
end
clc; %close all
addpath(genpath('MMEst'))
rmpath(genpath('MMEst/tucwav'))



X_orig{1} = cast(X_orig{1}, 'single');
X_orig{2} = cast(X_orig{2}, 'single');

[Y, ind_o, sigma] = genNm(X_orig{1}, rng_s);




% initializing parameters
[Rvph, Rtvs, Rtuc, ~] = initPart(tnsSize(Y), sigma, class(X_orig{1}));




% main loop
isTuc = true;
tic
[Pcs, ~, ~, nll] = Access(Y, ind_o, Rvph, Rtvs, Rtuc, isTuc, X_orig, []);
tim = toc;



X_hat = tnsMult(Pcs.S, Pcs.U_1, Pcs.U_2, Pcs.V);
X_hat3 = ROGCB2sRGB(X_hat);


rmpath(genpath('MMEst')); clc

end

