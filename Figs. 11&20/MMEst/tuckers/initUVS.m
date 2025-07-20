function [U_2, V, S] = initUVS(Y, isUV, DK)
%INITUVS Summary of this function goes here
%   Detailed explanation goes here


U_1 = [];
U_2 = [];
V = [];

if (isUV(1))
    Y_mde = tns2Mtx(Y, [1, 0, 0]);
    [U_1, ~, ~] = svd(Y_mde*Y_mde');
    U_1 = U_1(:, 1:DK(1));
end

if (isUV(2))
    Y_mde = tns2Mtx(Y, [0, 1, 0]);
    [U_2, ~, ~] = svd(Y_mde*Y_mde');
    U_2 = U_2(:, 1:DK(2));
end

if (isUV(3))
    Y_mde = tns2Mtx(Y, [0, 0, 1]);
    [V, ~, ~] = svd(Y_mde*Y_mde');
    %V = dct(eye(size(Y, 3), 'like', Y))';
    V = V(:, 1:DK(3));
end


S = tnsMult(Y, U_1', U_2', V');

end


