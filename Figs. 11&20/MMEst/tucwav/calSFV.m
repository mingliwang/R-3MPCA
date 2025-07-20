function nll = calSFV(theta, Y_wt, X_hat, mu, Fs, Etas, Rwv, S, isJTV, mg, Fg, Etag)
%CALSFV Summary of this function goes here
%   Detailed explanation goes here



Err = Y_wt - X_hat;
nll_1 = sum((Err(:)).^2);


SIZ = tnsSize(X_hat);
mdeMtx = logical(eye(3));
Xi_1 = calDiff(X_hat, SIZ, Fs.F_1, mdeMtx(1, :));
Xi_2 = calDiff(X_hat, SIZ, Fs.F_2, mdeMtx(2, :));
Xi_3 = calDiff(X_hat, SIZ, Fs.F_3, mdeMtx(3, :));
Res_1 = Xi_1 - Etas.Eta_1;
Res_2 = Xi_2 - Etas.Eta_2;
Res_3 = Xi_3 - Etas.Eta_3;
nll_2 = sum((Res_1(:)).^2) + sum((Res_2(:)).^2) + sum((Res_3(:)).^2);

if (isJTV)
    Xi_1 = 0;
    Xi_2 = 0;
    Xi_3 = calDiff(X_hat, SIZ, Fg.F_3, mdeMtx(3, :));
    Res_1 = Xi_1 - Etag.Eta_1;
    Res_2 = Xi_2 - Etag.Eta_2;
    Res_3 = Xi_3 - Etag.Eta_3;
    nll_2g = sum((Res_1(:)).^2) + sum((Res_2(:)).^2) + sum((Res_3(:)).^2);
else
    nll_2g = 0;
end


nll_3 = wvCoes(S, Rwv);


nll = (0.5*theta)*nll_1 + (0.5*mu)*nll_2 + nll_3 + (0.5*mg)*nll_2g;

end

