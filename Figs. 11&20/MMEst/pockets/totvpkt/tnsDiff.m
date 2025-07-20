function Xis = tnsDiff(X_hat, Fs)
%TNSDIFF Summary of this function goes here
%   Detailed explanation goes here



SIZ = tnsSize(X_hat);
mdeMtx = logical(eye(3));
Xis.Xi_1 = calDiff(X_hat, SIZ, Fs.F_1, mdeMtx(1, :));
Xis.Xi_2 = calDiff(X_hat, SIZ, Fs.F_2, mdeMtx(2, :));
Xis.Xi_3 = calDiff(X_hat, SIZ, Fs.F_3, mdeMtx(3, :));


end

