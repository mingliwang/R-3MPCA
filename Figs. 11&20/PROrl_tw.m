function X_hat = PROrl_tw(X_orig, rng_s)
%
%

for t = 1: 30
    maxNumCompThreads(1);
end
clc; %close all
addpath(genpath('MMEst'))
rmpath(genpath('MMEst/tuckers'))



X_orig = double(X_orig)/255;


[Y, ind_o, sigma] = genNm(X_orig, rng_s);
[Y, ind_o, SEPs] = tnsStre(Y, ind_o, 4);        % Rwv.scals = 4



% initializing parameters
[Rvph, Rtvs, ~, Rwv] = initPars(tnsSize(Y), sigma, class(X_orig));



% main loop
isTuc = false;
tic
[Pcs, ~, ~, nll] = Access(Y, ind_o, Rvph, Rtvs, Rwv, isTuc, X_orig, SEPs);
tim = toc;


X_hat = tnsMult(wvSyns(Pcs.S, Rwv), [], [], Pcs.V);
X_hat = tnsTrm(X_hat, SEPs);


rmpath(genpath('MMEst')); clc

end

