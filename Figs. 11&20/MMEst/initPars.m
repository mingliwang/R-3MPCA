function [Rvph, Rtvs, Rtuc, Rwv] = initPars(SIZ, sigma, varTyp)
%INITPARS Initializing Rvph, Rpsi, Rtuc
%   Detailed explanation goes here



% adjustable parameters
r_dp = 4e-2;
xis_dp = 9e-1;
gm_tv = 0.9e-2;

gm_s = 1e-32;                                                    % tucker, S_hat 1, 2
DK = [46, 46, 0];                                                 % for single input face
if (SIZ(3)>1)
    DK = [64, 64, 0];                                             % D_1, D_2, K; if 0, then bypass pc
end

if (strcmpi(varTyp, 'double'))
    gm_tv = 1.2e-2;                                              % for wavelet
end
gm_w = 1.5e-1;                                                  % wavelet, S_hat 1, 2
K = 0;


% varphi
gsm = 'h_4';
isBTF = false;                                                     % bundled tube-fiber
Rvph = initVph(cast(r_dp, varTyp), gsm, isBTF, SIZ(3));


% psi, F; psj
gsm = 'h_3';
isVTV = false;                                                     % vector TV
isMd3 = false;                                                     % psi, 1.1>true
isJTV = false;                                                      % psj, joined TV, for msi
Rtvs = initTVs(SIZ, cast(xis_dp, varTyp), gsm, isVTV, isMd3, gm_tv, sigma, isJTV);


% tucker
Rtuc = initTuc(DK, cast(gm_s, varTyp), sigma);


% wavelet
Rwv = initWv(K, gm_w, sigma);

end

