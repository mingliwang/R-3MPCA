function [V, S] = initVS(Y, isV, K, Rwv)
%INITUVS Summary of this function goes here
%   Detailed explanation goes here


V = [];

if (isV)
    Y_mde = tns2Mtx(Y, [0, 0, 1]);
    [V, ~, ~] = svd(Y_mde*Y_mde');
    %V = dct(eye(size(Y, 3)))';
    V = V(:, 1:K);
end


Z = tnsMult(Y, [], [], V');
S = wvAnas(Z, Rwv);

end


