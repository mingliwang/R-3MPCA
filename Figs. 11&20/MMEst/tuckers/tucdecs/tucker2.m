function Pcs = tucker2(Y, Rtuc, is0th, Pcs)
%TUCKER2 Summary of this function goes here
%   Detailed explanation goes here



% pre-variables
if (is0th)
    [U_2, S] = initUS(Y, Rtuc.DK);
else
    U_2 = Pcs.U_2;
    S = Pcs.S;
end

SIZ_Y = tnsSize(Y);
SIZ_S = [Rtuc.DK(1:2), SIZ_Y(3)];

Yt = permute(Y, [2, 1, 3]);


% 2d tucker decomposition
ITER = 2e0;
prec = 1e-6;

for iter = 1: ITER
    St = permute(S, [2, 1, 3]);
    U_1 = tdpca(Y, U_2, St, SIZ_Y(3));
    
    U_2 = tdpca(Yt, U_1, S, SIZ_Y(3));
    
    S_upd = tnsProd(Y, U_1, U_2, SIZ_S);
    if (mean(abs(S_upd(:)-S(:)))<prec)
        break
    end
    S = S_upd;
    
end

Pcs.U_1 = U_1;
Pcs.U_2 = U_2;
Pcs.V = [];
Pcs.S = S_upd;

end



function [U_2, S] = initUS(Y, DK)
%INITUS Summary of this function goes here
%   Detailed explanation goes here



Y_tmp = tns2Mtx(Y, [1, 0, 0]);
[U_1, ~, ~] = svd(Y_tmp*Y_tmp');
U_1 = U_1(:, 1:DK(1));

Y_tmp = tns2Mtx(Y, [0, 1, 0]);
[U_2, ~, ~] = svd(Y_tmp*Y_tmp');
U_2 = U_2(:, 1:DK(2));


S = tnsMult(Y, U_1', U_2', []);

end


function U_1 = tdpca(Y, U_2, St, N)
%TDPCA Summary of this function goes here
%   Detailed explanation goes here


B_1 = 0;
for n = 1: N
    Y_n = Y(:, :, n);
    St_n = St(:, :, n);
    
    B_1 = B_1 + (Y_n*U_2)*St_n;
end

[U_p, ~, V_p] = svd(B_1, 'econ');
U_1 = U_p*V_p';

end


function S = tnsProd(Y, U_1, U_2, SIZ_S)
%TNSPROD Summary of this function goes here
%   Detailed explanation goes here


U_1t = U_1';
S = zeros(SIZ_S, 'like', Y);
for n = 1: SIZ_S(3)
    S(:, :, n) = U_1t*Y(:, :, n)*U_2;
end

end


