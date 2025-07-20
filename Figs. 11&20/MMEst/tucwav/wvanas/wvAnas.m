function S = wvAnas(Z, Rwv)
%WVANAS Summary of this function goes here
%   Detailed explanation goes here



K = size(Z, 3);
S = cell(1, K);
scals = Rwv.scals;
Faf = Rwv.Faf;
af = Rwv.af;
for k = 1: K
    S{k} = cplxdual2D(Z(:, :, k), scals, Faf, af);
end

end

