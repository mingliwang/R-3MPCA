function Pcs = P_step(Y, Rvph, Rtvs, Rwv, Pcs, Tau, Oms)
%P_STEP Summary of this function goes here
%   Detailed explanation goes here



% extraction
theta = Rvph.theta;

Fs = Rtvs.Rpsi.Fs;
Ls = Rtvs.Rpsi.Ls;
varho = Rtvs.Rpsi.varho;
mu = Rtvs.Rpsi.mu;

isJTV = (~isempty(Rtvs.Rpsj));                         % joined TV
if (isJTV)
    Fg = Rtvs.Rpsj.Fs;
    Lg = Rtvs.Rpsj.Ls;
    varhg = Rtvs.Rpsj.varho;
    mg = Rtvs.Rpsj.mu;
else
    Fg = [];
    varhg = 0*varho;
    mg = 0*mu;
end


% pre-variables
Tau_nz = Tau/(Rvph.tau_1);

X_old = tnsMult(wvSyns(Pcs.S, Rwv), [], [], Pcs.V);
Y_wt = Tau_nz.*Y + (1-Tau_nz).*X_old;

Ome_nz = (Oms.Ome)/(Rtvs.Rpsi.ome_1);
[Etas, Ztas] = calZtas(X_old, Fs, Ome_nz);

if (isJTV)
    Omg_nz = (Oms.Omg)/(Rtvs.Rpsj.ome_1);
    [Etag, Ztag] = calZtas(X_old, Fg, Omg_nz);
else
    Etag = [];
    Ztag = [];
end

theta_tv = theta + mu*sum(varho) + mg*sum(varhg);
theta_nz = theta/theta_tv;
mu_nz = mu/theta_tv;
mg_nz = mg/theta_tv;


% iteration
ITER = 2;
prec = 1e-3;
nll_p = zeros(1, ITER, 'single');
nll_p(1) = calSFV(theta, Y_wt, X_old, mu, Fs, Etas, Rwv, Pcs.S, isJTV, mg, Fg, Etag);

X_oldp = X_old;
for iter = 1: ITER
    Et_ws = calEtws(X_oldp, Ls, varho, Ztas);
    if (isJTV)
        Et_wg = calEtws(X_oldp, Lg, varhg, Ztag);
    else
        Et_wg = 0;
    end
    
    Y_wtv = theta_nz*Y_wt + mu_nz*Et_ws + mg_nz*Et_wg;
    
    Pcs = tucker1(Y_wtv, Rwv, false, Pcs);
    Pcs.S = threshd(theta_tv, Rwv, Pcs);
    
    
    X_hat = tnsMult(wvSyns(Pcs.S, Rwv), [], [], Pcs.V);
    
    nll_p(iter+1) = calSFV(theta, Y_wt, X_hat, mu, Fs, Etas, Rwv, Pcs.S, isJTV, mg, Fg, Etag);
    if (abs(nll_p(iter+1)-nll_p(iter)) <= prec)
        break
    end
    X_oldp = X_hat;
    
end


end



