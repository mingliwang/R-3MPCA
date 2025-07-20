function S = threshd(theta_tv, Rtuc, Pcs)
%THRESHD Summary of this function goes here
%   Detailed explanation goes here



thrVal = (Rtuc.regl)/theta_tv;

S = Pcs.S;
if (Rtuc.isSof)
    sgn = sign(S);
    S = abs(S) - thrVal;
    S(S<0) = 0;
    S = S.*sgn;
else
    S(abs(S)<sqrt(2*thrVal)) = 0;
end


end

