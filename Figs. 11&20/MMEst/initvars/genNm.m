function [Y, ind_o, sigma] = genNm(X_orig, rng_s)
%GENNM Summary of this function goes here
%   Detailed explanation goes here



% adjustable parameters
per_obs = 1;                 % percentage of observation
per_spe = 0.6;                 % percentage of salt & pepper



% each frontal slice
N = size(X_orig, 3);
Y = 0*X_orig;
ind_o = 0*X_orig;
for n = 1: N
    X_n = X_orig(:, :, n);
    rng_t = rng_s + (n-1);
    
    [Y_n, ind_n] = genDm(per_obs, per_spe, X_n, rng_t);
    Y(:, :, n) = Y_n;
    ind_o(:, :, n) = ind_n;
end



% sigma
err = (Y - X_orig).*ind_o;
err = err(:);
sigma = sqrt(sum(err.^2)/sum(ind_o(:)));


% random initial
if (per_obs<1)
    Y = rdmInit(Y, ind_o);
end

end



function [Y, ind_o] = genDm(per_obs, per_spe, X, rng_t)
%GENDM Summary of this function goes here
%   Detailed explanation goes here



% pre-variable
[M_1, M_2] = size(X);
N_p = prod([M_1, M_2]);


% location
rng(rng_t)
loct = randperm(N_p);

N_o = floor(N_p*per_obs);                       % observed
N_os = floor((N_o*per_spe)/2);               % salt

loct_os = loct(1:N_os);                            % salt
loct_op = loct((N_os+1):(2*N_os));         % pepper


% ind of observation
ind_o = zeros(N_p, 1, 'like', X);
ind_o(loct(1:N_o)) = 1;
ind_o = reshape(ind_o, M_1, M_2);



% Y
Y = X(:);
Y(loct_os) = 1;
Y(loct_op) = 0;


Y = reshape(Y, M_1, M_2);
Y = Y.*ind_o;


end



function Y = rdmInit(Y, ind_o)
%RDMINIT Summary of this function goes here
%   Detailed explanation goes here



[M_1, M_2, N] = size(Y);


pval = [0, 60, 115, 170, 255]/255;
perc = [1, 2, 4, 2, 2];
perc = perc/sum(perc);


K = length(perc);
M = M_1*M_2;
num = floor(M*perc);
num(K) = M - sum(num(1:(K-1)));
RT = [];
for k = 1: K
    tmp = repmat(pval(k), 1, num(k));
    RT = cat(2, RT, tmp);
end


rng(0)
loct = randperm(M);
RT = RT(loct);
RT = reshape(RT, M_1, M_2);
RT = cast(RT, 'like', Y);
RT = repmat(RT, 1, 1, N);


Y = Y.*ind_o + RT.*(1-ind_o);

end



