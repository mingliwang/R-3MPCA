function [nll, Tau, Oms] = calOFV(Y, ind_o, Rvph, Rtvs, Rwv, Pcs)
%CALOFV Summary of this function goes here
%   Detailed explanation goes here



% X_hat; isJTV
X_hat = tnsMult(wvSyns(Pcs.S, Rwv), [], [], Pcs.V);
isJTV = (~isempty(Rtvs.Rpsj));                         % joined TV, for msi


% r_obs; tau, nll_phi
[Tau, nll_phi] = calNLL(Y, X_hat, ind_o, Rvph);


% ome, nll_psi
[Oms.Ome, nll_psi] = calTVs(X_hat, Rtvs.isVTV, Rtvs.Rpsi);
if (isJTV)
    [Oms.Omg, nll_psj] = calTVs(X_hat, Rtvs.isVTV, Rtvs.Rpsj);
else
    Oms.Omg = [];
    nll_psj = 0;
end


% S; regs
nll_3 = wvCoes(Pcs.S, Rwv);


% nll
coef_phi = (Rvph.theta)/(Rvph.tau_1);
coef_psi = (Rtvs.Rpsi.mu)/(Rtvs.Rpsi.ome_1);
if (isempty(Rtvs.Rpsj))
    coef_psj = 0;
else
    coef_psj = (Rtvs.Rpsj.mu)/(Rtvs.Rpsj.ome_1);
end
nll = coef_phi*(nll_phi) + coef_psi*(nll_psi) + coef_psj*(nll_psj) + nll_3;

end



function [Tau, nll_phi] = calNLL(Y, X_hat, ind_o, Rvph)
%CALNLL Summary of this function goes here
%   Detailed explanation goes here


isBTF = ((isfield(Rvph, 'isBTF'))&&(Rvph.isBTF));
if (isBTF)
    r_obs = (Y-X_hat).*ind_o + (Rvph.replc)*(1-ind_o);
    r_obs = sqrt(sum(r_obs.^2, 3));
    ind_o = [];
else
    r_obs = abs((Y-X_hat).*ind_o);
end

phi = Rvph.phi;  % phi, phid, r_dp, tau_1, beta_phi
phid = Rvph.phid;
r_dp = Rvph.r_dp;
tau_1 = Rvph.tau_1;
beta_phi = Rvph.beta_phi;

[Tau, nll_phi] = calNLF(phi, phid, r_dp, tau_1, beta_phi, r_obs, ind_o);
if (isBTF)
    Tau = repmat(Tau, 1, 1, size(Y, 3));
end

end



function [Ome, nll_psi] = calTVs(X_hat, isVTV, Rpsi)
%CALTVS Summary of this function goes here
%   Detailed explanation goes here



psi = Rpsi.psi;  % psi, psid, xi_dp, ome_1, beta_psi
psid = Rpsi.psid;
xis_dp = Rpsi.xis_dp;
ome_1 = Rpsi.ome_1;
beta_psi = Rpsi.beta_psi;

Xis = tnsDiff(X_hat, Rpsi.Fs);
Xi_fnm = (Xis.Xi_1).^2 + (Xis.Xi_2).^2 + (Xis.Xi_3).^2;

if (isVTV)
    Xi_fnm = sum(Xi_fnm, 3);
end
Xi_fnm = sqrt(Xi_fnm);

[Ome, nll_psi] = calNLF(psi, psid, xis_dp, ome_1, beta_psi, Xi_fnm, []);
if (isVTV)
    Ome = repmat(Ome, 1,  1, size(X_hat, 3));
end

end



function [tau, nll_phi] = calNLF(phi, phid, r_dp, tau_1, beta_phi, r_obs, ind_o)
%CALNLF Summary of this function goes here
%   Detailed explanation goes here



if (isempty(ind_o))
    ind_o = ones(1, 'like', r_obs);
end

ind_1 = (abs(r_obs)<=r_dp).*ind_o;
ind_2 = (1 - ind_1).*ind_o;

s_obs = (abs(r_obs).^2)/2;

tau = tau_1*ind_1 + (phid(s_obs)).*ind_2;

nll_1 = tau_1*(s_obs);
nll_2 = phi(s_obs) + beta_phi;
nll_phi = nll_1.*ind_1 + nll_2.*ind_2;
nll_phi = sum(nll_phi(:));


end


