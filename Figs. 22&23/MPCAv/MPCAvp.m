function [uptInRom, prsInRom] = MPCAvp(stdGlb, ptsInRom, pteInRom, simIndNbh, patSize, simNum, balChl, balWie)
%MPCAVP Denoising function
%   Y(:, column, :) is a patch


uptInRom = zeros(size(ptsInRom), 'single');
prsInRom = uptInRom;

balGN = single(0.9);
hrdThr = single(0.5);  % hard-threshold for shrinkage coefficient
for krc = 1: size(simIndNbh, 1)
    curInds = simIndNbh(krc, :);
    Y = ptsInRom(:, curInds);
    Z = pteInRom(:, curInds);
    
    % remove mean, and get varGN&prsGN
    Ym = repmat(sum(Y, 2)/simNum, 1, simNum);
    Yr = Y - Ym;
    stdNbh = sqrt(sum(sum(Yr.^2, 1), 2)/((patSize(1)*patSize(1)*patSize(2)*simNum)-1));
    varGN = (balGN*stdGlb+(1-balGN)*stdNbh).^2;
    prsGN = repmat((((1-balGN)*stdGlb+balGN*stdNbh).^(-1)).^2, patSize(1)*patSize(1)*patSize(2), simNum);
    
    % pilot for denoising
    Zm = repmat(sum(Z, 2)/simNum, 1, simNum);
    Zr = Z - Zm;
    [U, ~] = eig(balChl*(Zr*Zr') + (Yr*Yr'));
    ZtaS = sum((U'*Zr).^2, 2);
    
    rShc = ZtaS.*((ZtaS+balWie*varGN).^(-1));
    rShc = single(rShc>=hrdThr);
    X = (U*diag(rShc)*U')*Yr + Ym;
    
    uptInRom(:, curInds) = uptInRom(:, curInds) + X.*prsGN;
    prsInRom(:, curInds) = prsInRom(:, curInds) + prsGN;
end
end

