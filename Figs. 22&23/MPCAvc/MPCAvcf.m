function [uptInRom, prsInRom] = MPCAvcf(stdGlb, ptsInRom, simIndNbh, patSize, simNum, balDiv)
%MPCAVCF Denoising function, for the first principal channel
%   Y(:, column, :) is a patch


uptInRom = zeros(size(ptsInRom), 'single');
prsInRom = uptInRom;

balGN = single(0.9);
hrdThr = single(0.5);  % hard-threshold for shrinkage coefficient
for krc = 1: size(simIndNbh, 1)
    curInds = simIndNbh(krc, :);
    Y = ptsInRom(:, curInds, :);
    
    % remove mean, and get varGN&prsGN
    Ym = repmat(sum(Y, 2)/simNum, 1, simNum);
    Yr = Y - Ym;
    stdNbh = sqrt(sum(sum(Yr.^2, 1), 2)/((patSize(1)*patSize(1)*patSize(2)*simNum)-1));
    varGN = (balGN*stdGlb+(1-balGN)*stdNbh).^2;
    prsGN = repmat((((1-balGN)*stdGlb+balGN*stdNbh).^(-1)).^2, patSize(1)*patSize(1)*patSize(2), simNum, 1);
    
    
    % for denoising
    X = Y;
    Yrf = Yr(:, :, 1);
    [U, ~] = eig(Yrf*Yrf');
    Eta = U'*Yrf;
    rShc = 1 - balDiv(1)*varGN(1)*((sum(Eta.^2, 2)).^(-1));
    rShc = single(rShc>=hrdThr);
    X(:, :, 1) = U*(repmat(rShc, 1, simNum).*Eta) + Ym(:, :, 1);
    
    uptInRom(:, curInds, :) = uptInRom(:, curInds, :) + X.*prsGN;
    prsInRom(:, curInds, :) = prsInRom(:, curInds, :) + prsGN;
end
end

