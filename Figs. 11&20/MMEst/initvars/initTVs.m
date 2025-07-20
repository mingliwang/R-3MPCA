function Rtvs = initTVs(SIZ, xis_dp, gsm, isVTV, isMd3, gm_tv, sigma, isJTV)
%INITTVS Summary of this function goes here
%   Detailed explanation goes here



order = 1;
varTyp = class(xis_dp);
Fs.F_1 = circMtx(SIZ(1), order, varTyp);
Fs.F_2 = circMtx(SIZ(2), order, varTyp);
Fs.F_3 = circMtx(SIZ(3), order, varTyp);
Ls.L_1 = (Fs.F_1)'*(Fs.F_1);
Ls.L_2 = (Fs.F_2)'*(Fs.F_2);
Ls.L_3 = (Fs.F_3)'*(Fs.F_3);
varho = (sum(abs(Fs.F_1(1, :))))^2;


Rpsi = initPsi(SIZ, xis_dp, gsm, isVTV, isMd3, gm_tv, sigma, Fs, Ls, varho);

if (isJTV)
    Rpsj = initPsj(SIZ, xis_dp, gsm, isVTV, gm_tv, sigma, Fs, Ls, varho);
else
    Rpsj = [];
end

Rtvs.isVTV = isVTV;
Rtvs.Rpsi = Rpsi;
Rtvs.Rpsj = Rpsj;

end



function Rpsi = initPsi(SIZ, xis_dp, gsm, isVTV, isMd3, gm_tv, sigma, Fs, Ls, varho)
%INITPSI Summary of this function goes here
%   Detailed explanation goes here



varTyp = class(xis_dp);

if (isVTV)
    N = SIZ(3);
    xis_dp = xis_dp*sqrt(N);
else
    N = cast(1, varTyp);
end
isMd3 = cast(isMd3, varTyp);
Np = N*((2+isMd3^2)/2);
[psi, psid] = gsmLst(gsm, Np);

s_dp = (xis_dp(1)^2)/2;
ome_1 = psid(s_dp);
beta_psi = ome_1*(s_dp) - psi(s_dp);


Rpsi.psi = psi;
Rpsi.psid = psid;
Rpsi.xis_dp = xis_dp(1);
Rpsi.ome_1 = ome_1;
Rpsi.beta_psi = beta_psi;
Rpsi.Fs = Fs;
Rpsi.Ls = Ls;
Rpsi.varho = varho*[1, 1, 1];


if (isMd3==0)
    Rpsi.Fs.F_3 = cast(0, varTyp);
    Rpsi.Ls.L_3 = cast(0, varTyp);
else
    Rpsi.Fs.F_3 = isMd3*(Rpsi.Fs.F_3);
    Rpsi.Ls.L_3 = (isMd3^2)*(Rpsi.Ls.L_3);
end
Rpsi.varho(3) = (isMd3^2)*(Rpsi.varho(3));

Rpsi.mu = cast(gm_tv(1)*((sqrt(112*92)/4)*sigma), varTyp);

end



function Rpsj = initPsj(SIZ, xis_dp, gsm, isVTV, gm_tv, sigma, Fs, Ls, varho)
%INITPSJ Summary of this function goes here
%   Detailed explanation goes here



varTyp = class(xis_dp);

if (length(xis_dp)==1)
    xis_dp(2) = cast(1e-1, varTyp);
end

if (isVTV)
    N = SIZ(3);
    xis_dp = xis_dp*sqrt(N);
else
    N = cast(1, varTyp);
end
[psi, psid] = gsmLst(gsm, N/2);
s_dp = (xis_dp(2)^2)/2;
ome_1 = psid(s_dp);
beta_psi = ome_1*(s_dp) - psi(s_dp);


Rpsj.psi = psi;
Rpsj.psid = psid;
Rpsj.xis_dp = xis_dp(2);
Rpsj.ome_1 = ome_1;
Rpsj.beta_psi = beta_psi;
Rpsj.Fs.F_1 = cast(0, varTyp);
Rpsj.Fs.F_2 = cast(0, varTyp);
Rpsj.Fs.F_3 = Fs.F_3;
Rpsj.Ls.L_1 = cast(0, varTyp);
Rpsj.Ls.L_2 = cast(0, varTyp);
Rpsj.Ls.L_3 = Ls.L_3;
Rpsj.varho = varho*[0, 0, 1];


if (length(gm_tv)==1)
    gm_tv(2) = gm_tv(1);
end
Rpsj.mu = cast(gm_tv(2)*((sqrt(112*92)/4)*sigma), varTyp);

end


