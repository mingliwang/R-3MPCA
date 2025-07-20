function S_hat = threshd(theta_tv, Rwv, Pcs)
%THRESHD Summary of this function goes here
%   Detailed explanation goes here



scals = Rwv.scals;
isSof = Rwv.isSof;

thrVal = (Rwv.regs)/theta_tv;


S = Pcs.S;
K = length(S);
S_hat = cell(1, K);
I = sqrt(-1);
for k = 1: K
    w = S{k};
    for j = 1: scals
        for s1 = 1: 2
            for s2 = 1: 3
                C = w{j}{1}{s1}{s2} + I*w{j}{2}{s1}{s2};
                if (isSof)
                    C = soft(C, thrVal(j));
                else
                    ind = (abs(C)>thrVal(j));
                    C = ind.*C;
                end
                w{j}{1}{s1}{s2} = real(C);
                w{j}{2}{s1}{s2} = imag(C);
            end
        end
    end
    S_hat{k} = w;
end


end



