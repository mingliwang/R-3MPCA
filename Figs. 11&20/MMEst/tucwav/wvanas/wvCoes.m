function coeVal = wvCoes(S, Rwv)
%WVCOES Summary of this function goes here
%   Detailed explanation goes here



scals = Rwv.scals;
isSof = Rwv.isSof;
regs = Rwv.regs;

K = length(S);
coeVal = zeros(K, scals);

I = sqrt(-1);
for k = 1: K
    w = S{k};
    
    for j = 1: scals
        tmp = 0;
        for s1 = 1:2
            for s2 = 1:3
                C = w{j}{1}{s1}{s2} + I*w{j}{2}{s1}{s2};
                if (isSof)
                    tmp = tmp + sum(abs(C(:)));
                else
                    tmp = tmp + sum(abs(C(:))>1e-12);
                end
            end
        end
        coeVal(k, j) = tmp;
        
    end
end


coeVal = coeVal.*repmat(regs, K, 1);
coeVal = sum(coeVal(:));

end
