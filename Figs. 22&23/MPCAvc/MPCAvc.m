function [uptInRom, prsInRom] = MPCAvc(stdGlb, ptsInRom, simIndNbh, patSize, simNum, balChl, balDiv, isBnr)
%MPCAVC Denoising function
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
    [U, ~] = eig(balChl(1)*(Yr(:, :, 1)*Yr(:, :, 1)')+balChl(2)*(Yr(:, :, 2)*Yr(:, :, 2)')+balChl(3)*(Yr(:, :, 3)*Yr(:, :, 3)'));
    for dim = 1: 3
        Eta = U'*Yr(:, :, dim);
        rShc = 1 - balDiv(dim)*varGN(dim)*((sum(Eta.^2, 2)).^(-1));
        if isBnr
            rShc = single(rShc>=hrdThr);
        else
            rShc = max(rShc, single(rShc>=hrdThr));
        end
        X(:, :, dim) = U*(repmat(rShc, 1, simNum).*Eta) + Ym(:, :, dim);
    end
    
    uptInRom(:, curInds, :) = uptInRom(:, curInds, :) + X.*prsGN;
    prsInRom(:, curInds, :) = prsInRom(:, curInds, :) + prsGN;
end
end

