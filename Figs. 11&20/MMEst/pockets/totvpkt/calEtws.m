function Ewts = calEtws(X_hat, Ls, varho, Ztas)
%CALETWS Summary of this function goes here
%   Detailed explanation goes here



SIZ = tnsSize(X_hat);
mdeMtx = logical(eye(3));
Eta_1wt = varho(1)*(X_hat) - calDiff(X_hat, SIZ, Ls.L_1, mdeMtx(1, :)) + Ztas.Zta_1;  % varho(m)*R_mwt
Eta_2wt = varho(2)*(X_hat) - calDiff(X_hat, SIZ, Ls.L_2, mdeMtx(2, :)) + Ztas.Zta_2;
Eta_3wt = varho(3)*(X_hat) - calDiff(X_hat, SIZ, Ls.L_3, mdeMtx(3, :)) + Ztas.Zta_3;


Ewts = Eta_1wt + Eta_2wt + Eta_3wt;


end

