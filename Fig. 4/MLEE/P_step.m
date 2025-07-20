function [U, nll_p] = P_step(Y, Plin, tau, U)
%P_STEP PCA step
%   MM for updating A.


% double type
tau = double(tau);

% pre-variables
K = Plin.K;
r_sw = Plin.r_sw;
r_e = Plin.r_e;
A_l = Plin.A_l;
A_u = Plin.A_u;

[M, N] = size(Y);
varsig = (Plin.A_ui)*tau;

Z_igs = calNmc(r_e, M, K, Plin.A_lui, tau);  % Z value


% update A
ITER = 1e3;
prec = 1e-6;
nll_p = zeros(1, ITER, 'double');
for iter = 1: ITER
    r = double(sqrt(sum((Y-(U*U')*Y).^2, 1)));
    rLoct = double(calLoct(r, r_sw, K));
    
    nll_p(iter) = calNLL(Z_igs, r, rLoct, N, A_l, A_u, varsig);
    if ((iter>=2)&&(abs(nll_p(iter)-nll_p(iter-1))<prec))
        break
    end
    
    [U, ~, ~] = svds(Y*diag(tau(rLoct))*Y', 1);
end

nll_p = nll_p(1:iter);
end

