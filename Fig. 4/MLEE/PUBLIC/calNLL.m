function nll = calNLL(Z_igs, r, rLoct, N, A_l, A_u, varsig)
%CALNLL Summary of this function goes here
%   Detailed explanation goes here


% value of sum_{n, k} h_k(r_n^2/2)
h_val = r.*0;
for n = 1: N
    tmp = A_l(rLoct(n), :) + A_u(rLoct(n), :)*((r(n).^2)/2);
    h_val(n) = tmp*varsig;
end
h_val = sum(h_val);


% nll value, omit -sum_{n} log(r_n^(M-1))
nll = N*log(Z_igs) + h_val;  % F value, Psi=I, omit N*ln(2*(pi^{M/2})/Gam(M/2))
end

