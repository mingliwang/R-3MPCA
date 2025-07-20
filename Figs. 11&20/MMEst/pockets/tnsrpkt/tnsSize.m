function SIZ = tnsSize(X)
%TNSSIZE Summary of this function goes here
%   Detailed explanation goes here



[M_1, M_2, N] = size(X);
SIZ = cast([M_1, M_2, N], 'like', X);


end

