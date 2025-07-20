function [Rvph, Rtvs, Rtuc, Rwv] = initPart(SIZ, sigma, varTyp)
%INITPART Initializing Rvph, Rpsi, Rtuc
%   Detailed explanation goes here



% adjustable parameters
r_dp = 2.6e-2;
xis_dp = 2e-1;
gm_tv = 0.6e-2;


gm_s = 1e-32;                                                    % tucker, S_hat 1, 2
DK = floor([40, 40, 3]);                                       % 56, 40; 3; if 0, then bypass pc

gm_w = 3e-2;                                                     % wavelet, S_hat 1, 2
K = 0;



% varphi
gsm = 'h_4';
isBTF = false;                                                     % bundled tube-fiber
Rvph = initVph(cast(r_dp, varTyp), gsm, isBTF, SIZ(3));


% psi, F; psj
gsm = 'h_3';
isVTV = true;                                                     % vector TV
isMd3 = false;                                                     % psi, 1.1>true
isJTV = false;                                                      % psj, joined TV, for msi
Rtvs = initTVs(SIZ, cast(xis_dp, varTyp), gsm, isVTV, isMd3, gm_tv, sigma, isJTV);


% tucker
Rtuc = initTuc(DK, cast(gm_s, varTyp), sigma);


% wavelet
Rwv = initWv(K, gm_w, sigma);

end

