function Z = wvSyns(S, Rwv)
%WVSYNS Summary of this function goes here
%   Detailed explanation goes here



scals = Rwv.scals;
Fsf = Rwv.Fsf;
sf = Rwv.sf;
K = length(S);


Z = icplxdual2D(S{1}, scals, Fsf, sf);
Z = repmat(Z, 1, 1, K);
if (K>=2)
    for k = 2: K
        Z(:, :, k) = icplxdual2D(S{k}, scals, Fsf, sf);
    end
end

end
