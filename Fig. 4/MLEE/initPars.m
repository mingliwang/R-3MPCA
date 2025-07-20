function [U, Plin, tau] = initPars(Y)
%INITPARS Summary of this function goes here
%   Detailed explanation goes here


[U, ~, ~] = svds(Y*Y', 1);

K = 9;  % tau_{1:K}


% r_e
N_K = 4;  % number of points in (s_{K-1}, 1e256)
err = Y-(U*U')*Y;
[r_sw, r_e] = er2Edgs(err, K, N_K);                  % this style can be modified


% tau
tau = zeros(K, 1, 'double');
if (K==1)
    tau = 1;
else
    tau(1:(K-1)) = ((r_e(2:K)+r_e(1:(K-1)))/2).^(-2);  % t-type
    tau(K) = tau(K-1)/8;
end



% h = (A_l+A_u*s)*varsig
A_l = zeros(K, K, 1, 'double');
for k = 1: (K-1)
    A_l((k+1):K, k) = (r_e(k+1).^2)/2;
end
A_u = triu(ones(K, K, 'double'));  % tau=A_u*varsig
A_ui = A_u\eye(K);  % varsig=A_ui*tau
A_lui = A_l*A_ui;  % for calNmc

% Plin structure
Plin.K = K;
Plin.r_sw = double(r_sw);
Plin.r_e = double(r_e);
Plin.A_l = double(A_l);
Plin.A_u = double(A_u);
Plin.A_ui = double(A_ui);
Plin.A_lui = double(A_lui);
end


function [r_sw, r_e] = er2Edgs(err, K, N_K)
%ER2EDGS Partition points of R^+
%   Bin, histogram.


r = sqrt(sum(err.^2, 1));
N = length(r);
if (K>N)
    K = N;
end

% s_0=0, s_1, ..., s_K=1e256
r_sort = sort(r, 'ascend');
r_eK = (r_sort(N-N_K)+r_sort(N-N_K+1))/2;
r_sw = r_eK/(K-1);  % step width
r_e = (1:1:(K-1))*r_sw;  % edge
if (N==1)
    r_e = [0, 1e256];
else
    r_e = cat(2, 0, r_e, 1e256);
end
end

