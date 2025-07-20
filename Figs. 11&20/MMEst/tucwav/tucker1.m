function Pcs = tucker1(Y, Rwv, is0th, Pcs)
%TUCKER1 Summary of this function goes here
%   Detailed explanation goes here


% pre-variables
isV = Rwv.isV;
K = Rwv.K;

if (is0th)
    [V, S] = initVS(Y, isV, K, Rwv);
else
    V = Pcs.V;
    S = Pcs.S;
end


% tuc1+wav2 decomposition
ITER = 1;

for iter = 1: ITER
    V = pca1d(Y, isV, Rwv, S);
    
    S = wav2d(Y, isV, Rwv, V);
    
end


Pcs.V = V;
Pcs.S = S;

end



function V = pca1d(Y, isV, Rwv, S)
%PCA1D Summary of this function goes here
%   Detailed explanation goes here


if (~isV)
    V = [];
    return
end


Z = wvSyns(S, Rwv);


SIZ = tnsSize(Y);
Y = reshape(Y, SIZ(1)*SIZ(2), SIZ(3));
Y = Y';
Z = reshape(Z, SIZ(1)*SIZ(2), length(S));

B = Y*Z;
[U_B, ~, V_B] = svd(B, 'econ');
V = U_B*V_B';


end



function S = wav2d(Y, isV, Rwv, V)
%WAV2D Summary of this function goes here
%   Detailed explanation goes here



if (isV)
    Z = tnsMult(Y, [], [], V');
else
    Z = Y;
end


S = wvAnas(Z, Rwv);

end


