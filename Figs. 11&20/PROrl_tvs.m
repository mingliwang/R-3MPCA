function X_hat = PROrl_tvs(X_orig, rng_s)
%
%

for t = 1: 30
    maxNumCompThreads(1);
end
clc; %close all
addpath(genpath('MMEst'))
rmpath(genpath('MMEst/tucwav'))



X_orig = single(X_orig)/255;


[Y, ind_o, sigma] = genNm(X_orig, rng_s);




% initializing parameters
[Rvph, Rtvs, Rtuc, ~] = initPars(tnsSize(Y), sigma, class(X_orig));



% main loop
isTuc = true;
tic
[Pcs, ~, ~, nll] = Access(Y, ind_o, Rvph, Rtvs, Rtuc, isTuc, X_orig, []);
tim = toc;


X_hat = tnsMult(Pcs.S, Pcs.U_1, Pcs.U_2, Pcs.V);



rmpath(genpath('MMEst')); clc

end

