function Rwv = initWv(K, gm_w, sigma)
%INITWV Summary of this function goes here
%   Detailed explanation goes here



varTyp = 'double';

Rwv.isV = logical(K);
Rwv.K = cast(K, varTyp);


Rwv.scals = 4;  % number of scales

[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;

Rwv.Faf = Faf;
Rwv.Fsf = Fsf;
Rwv.af = af;
Rwv.sf = sf;


Rwv.isSof = false;
Rwv.regs = cast(gm_w*(sigma^2)*(1./(1.^(-1+(1:1:(Rwv.scals))))), varTyp);       % 1^()--2^()

end

