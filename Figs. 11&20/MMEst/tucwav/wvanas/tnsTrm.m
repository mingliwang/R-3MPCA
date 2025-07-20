function Y = tnsTrm(Y_etr, SEPs)
%TNSTRM Summary of this function goes here
%   Detailed explanation goes here




Y = Y_etr(SEPs(1, 1, 1):SEPs(1, 2, 1), ...
    SEPs(1, 1, 2):SEPs(1, 2, 2), :);



end

