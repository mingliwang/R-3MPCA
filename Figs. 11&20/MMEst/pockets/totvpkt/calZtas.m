function [Etas, Ztas] = calZtas(X_old, Fs, Ome_nz)
%CALZTAS Summary of this function goes here
%   Detailed explanation goes here



SIZ = tnsSize(X_old);
mdeMtx = logical(eye(3));
Etas.Eta_1 = (1-Ome_nz).*(calDiff(X_old, SIZ, Fs.F_1, mdeMtx(1, :)));    % R_m
Etas.Eta_2 = (1-Ome_nz).*(calDiff(X_old, SIZ, Fs.F_2, mdeMtx(2, :)));
Etas.Eta_3 = (1-Ome_nz).*(calDiff(X_old, SIZ, Fs.F_3, mdeMtx(3, :)));

Ztas.Zta_1 = calDiff(Etas.Eta_1, SIZ, (Fs.F_1)', mdeMtx(1, :));                  % R_m *_m F_m'
Ztas.Zta_2 = calDiff(Etas.Eta_2, SIZ, (Fs.F_2)', mdeMtx(2, :));
Ztas.Zta_3 = calDiff(Etas.Eta_3, SIZ, (Fs.F_3)', mdeMtx(3, :));


end

