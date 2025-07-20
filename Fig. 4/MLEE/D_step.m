function [tau, nll_d] = D_step(Y, Plin, U, tau)
%D_STEP Density estimation step
%   Max likelihood estimate for tau_{1:K}.



% double type
tau = double(tau);


% pre-variables
K = Plin.K;
r_e = Plin.r_e;
A_l = Plin.A_l;
A_u = Plin.A_u;
A_ui = Plin.A_ui;
A_lui = Plin.A_lui;

[M, N] = size(Y);
r = double(sqrt(sum((Y-(U*U')*Y).^2, 1)));
rLoct = double(calLoct(r, Plin.r_sw, K));


% h_grad; Z_pt1, Z_pp1
h_grad = calHgd(r, rLoct, N, K, A_l, A_u);
Z_pt1 = -(A_l+(1/2)*A_u)';
Z_pp1 = (A_l.^2+(1/4)*A_u)';


% update varsigj (element-wise)
ITER = 1e3;
prec = 1e-6;
nll_d = zeros(1, ITER, 'double');
for iter = 1: ITER  % for tau
    Z_igs = calNmc(r_e, M, K, A_lui, tau);  % Z value
    
    nll_d(iter) = calNLL(Z_igs, r, rLoct, N, A_l, A_u, A_ui*tau);
    if ((iter>=2)&&(abs(nll_d(iter)-nll_d(iter-1))<prec))
        break
    end
    
    varsig = A_ui*tau;
    for j = 1:1:K  % j = K:(-1):1
        varsig(j) = NewRp(r_e, M, K, A_l, A_u, h_grad, Z_pt1(j, :), Z_pp1(j, :), prec, N/Z_igs, varsig, j);
    end
    tau = A_u*varsig;
end

nll_d = nll_d(1:iter);
end



function h_grad = calHgd(r, rLoct, N, K, A_l, A_u)
%CALHGD Calculate gradient of h
%   Detailed explanation goes here


A_lt = A_l';
A_ut = A_u';

h_grad = zeros(K, 1, 'double');
for n = 1: N
    tmp = A_lt(:, rLoct(n)) + A_ut(:, rLoct(n))*((r(n).^2)/2);  % rLoct(n) determines h_k
    h_grad = h_grad + tmp;
end
end


