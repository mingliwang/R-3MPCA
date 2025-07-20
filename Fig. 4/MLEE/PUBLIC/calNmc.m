function Z_igs = calNmc(r_e, M, K, A_lui, tau)
%CALNMC Calculate normalized coefficient
%   r_e(1)=0, r_e(end)=1e256



rho = A_lui*tau;  % sum_{1...k-1} varsigi*si

Z_igs = zeros(1, K, 'double');
for k = 1: K
    Z_igs(k) = exp(-rho(k))*gmsinc(r_e(k), r_e(k+1), M, tau(k));
end

Z_igs = sum(Z_igs);
end
