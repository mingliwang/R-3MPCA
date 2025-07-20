function P_pca = PROw(outNum)
%
%


% clear; clc; close all
for t = 1: 30
    maxNumCompThreads(1);
end

addpath(genpath('MLEE'))


% load data
Y = loadOF(outNum);

[U, Plin, tau] = initPars(Y);


P_pca = cell(2, 3);
P_pca{1, 1} = sprintf('outVals');
P_pca{2, 1} = Y(:, (end-outNum+1):end);
P_pca{1, 2} = 'Plin';
P_pca{2, 2} = Plin;
P_pca{1, 3} = 'U_rec';
P_pca{1, 4} = 'Tau_rec';
P_pca{1, 5} = 'Nll_rec';
P_pca{1, 6} = 'LOOP';
P_pca{1, 7} = 'Time(s)';
P_pca{1, 8} = 'Lop1nll_d';

% D-step, P_step
LOOP = 2e2;
prec = 1e-6;

Tau_rec = zeros(length(tau), LOOP+1, 'single');
Tau_rec(:, 1) = single(tau);
U_rec = zeros(size(U, 1), LOOP+1, 'single');
U_rec(:, 1) = single(U);
Nll_rec = zeros(2, LOOP, 'single');

tic
for loop = 1: LOOP
    % PCA step
    [U, nll_p] = P_step(Y, Plin, tau, U);
    %plot(nll_p)
    
    % density estimation step
    [tau, nll_d] = D_step(Y, Plin, U, tau);
    %plot(nll_d)
    
    if (loop==1)
        nll_0(1, 1) = single(0);
        nll_0(2, 1) = single(nll_p(1));
        
        P_pca{2, 8} = single(nll_d);
    end
    
    U_rec(:, loop+1) = single(U);
    Tau_rec(:, loop+1) = single(tau);
    Nll_rec(1, loop) = single(nll_p(end));
    Nll_rec(2, loop) = single(nll_d(end));
    
    if (abs(nll_p(end)-nll_d(end))<=prec)
        break
    end
end
tim = toc;

P_pca{2, 3} = single(U_rec(:, 1:(loop+1)));
P_pca{2, 4} = single(Tau_rec(:, 1:(loop+1)));
P_pca{2, 5} = single(cat(2, nll_0, Nll_rec(:, 1:(loop))));
P_pca{2, 6} = single(loop);
P_pca{2, 7} = single(tim);


rmpath(genpath('MLEE'))
end


