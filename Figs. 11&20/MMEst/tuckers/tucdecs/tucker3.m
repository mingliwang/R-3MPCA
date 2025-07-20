function Pcs = tucker3(Y, Rtuc, is0th, Pcs)
%TUCKER3 Summary of this function goes here
%   Detailed explanation goes here


% pre-variables
isUV = Rtuc.isUV;
SIZ_Y = tnsSize(Y);
SIZ_S = Rtuc.DK;

if (is0th)
    [U_2, V, S] = initUVS(Y, isUV, SIZ_S);
else
    U_2 = Pcs.U_2;
    V = Pcs.V;
    S = Pcs.S;
end


% 2d tucker decomposition
ITER = 2e0;
prec = 1e-6;

for iter = 1: ITER
    [U_1, U_2] = pca2d(Y, isUV, U_2, V, S, SIZ_Y, SIZ_S(3));
    V = pca1d(Y, isUV, U_1, U_2, S);
    
    S_upd = tnsMult(Y, U_1', U_2', V');
    if (mean(abs(S_upd(:)-S(:)))<prec)
        break
    end
    S = S_upd;
    
end

Pcs.U_1 = U_1;
Pcs.U_2 = U_2;
Pcs.V = V;
Pcs.S = S_upd;

end



function [U_1, U_2] = pca2d(Y, isUV, U_2, V, S, SIZ, K)
%PCA2D Summary of this function goes here
%   Detailed explanation goes here


if (isUV(3))
    Y = reshape(Y, SIZ(1)*SIZ(2), SIZ(3));
    SIZ(3) = K;
    Y = reshape(Y*V, SIZ);
end

if (isUV(1))
    B = 0;
    St = permute(S, [2, 1, 3]);
    for k = 1: SIZ(3)
        Y_k = Y(:, :, k);
        St_k = St(:, :, k);
        if (isUV(2))
            B = B + (Y_k*U_2)*St_k;
        else
            B = B + Y_k*St_k;
        end
    end
    [U_B, ~, V_B] = svd(B, 'econ');
    U_1 = U_B*V_B';
else
    U_1 = [];
end

if (isUV(2))
    B = 0;
    Y = permute(Y, [2, 1, 3]);
    for k = 1: SIZ(3)
        Y_k = Y(:, :, k);
        S_k = S(:, :, k);
        if (isUV(1))
            B = B + (Y_k*U_1)*S_k;
        else
            B = B + Y_k*S_k;
        end
    end
    [U_B, ~, V_B] = svd(B, 'econ');
    U_2 = U_B*V_B';
end

end



function V = pca1d(Y, isUV, U_1, U_2, S)
%PCA1D Summary of this function goes here
%   Detailed explanation goes here


if (~isUV(3))
    V = [];
    return
end


Y = tnsMult(Y, U_1', U_2', []);

SIZ = tnsSize(Y);
Y = reshape(Y, SIZ(1)*SIZ(2), SIZ(3));
Y = Y';

S = reshape(S, SIZ(1)*SIZ(2), []);

B = Y*S;
[U_B, ~, V_B] = svd(B, 'econ');
V = U_B*V_B';


end


