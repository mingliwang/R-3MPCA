function [uptInRom, prsInRom] = MPCAv(stdGlb, ptsInRom, simIndNbh, patSize, simNum, balDiv, isBnr)
%MPCAV Denoising function
%   Y(:, column, :) is a patch


uptInRom = zeros(size(ptsInRom), 'single');
prsInRom = uptInRom;

balGN = single(0.9);
hrdThr = single(0.5);  % hard-threshold for shrinkage coefficient
for krc = 1: size(simIndNbh, 1)
    curInds = simIndNbh(krc, :);
    Y = ptsInRom(:, curInds);
    
    % remove mean, and get varGN&prsGN
    Ym = repmat(sum(Y, 2)/simNum, 1, simNum);
    Yr = Y - Ym;
    stdNbh = sqrt(sum(sum(Yr.^2, 1), 2)/((patSize(1)*patSize(1)*patSize(2)*simNum)-1));
    varGN = (balGN*stdGlb+(1-balGN)*stdNbh).^2;
    prsGN = repmat((((1-balGN)*stdGlb+balGN*stdNbh).^(-1)).^2, patSize(1)*patSize(1)*patSize(2), simNum);
    
    % for denoising
    [U, S2] = eig(Yr*Yr');
    sShc = 1 - (balDiv*varGN)*(abs(diag(S2)).^(-1));
    if isBnr
        sShc = single(sShc>=hrdThr);
    else
        sShc = max(sShc, single(sShc>=hrdThr));
    end
    X = (U*diag(sShc)*U')*Yr + Ym;
    
    uptInRom(:, curInds) = uptInRom(:, curInds) + X.*prsGN;
    prsInRom(:, curInds) = prsInRom(:, curInds) + prsGN;
end
end

