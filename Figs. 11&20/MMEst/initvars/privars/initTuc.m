function Rtuc = initTuc(DK, gm_s, sigma)
%INITTUC Summary of this function goes here
%   Detailed explanation goes here



varTyp = class(gm_s);

Rtuc.isUV = logical(DK);
Rtuc.DK = cast(DK, varTyp);
Rtuc.isSof = true;
Rtuc.regl = cast(gm_s*(sigma^2), varTyp);

end

