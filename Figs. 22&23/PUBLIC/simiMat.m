function simIndNbh = simiMat(stdGlb, simNum, ptsInRom, rptIndNbh)
%SIMIMAT Matching similar patches for reference patches
%   by using Euclidean norm.

if (size(ptsInRom, 3)>1)
    ptsInRom = ptsInRom(:, :, 1);  % using first principal channel
end

% denoising step
if (stdGlb>30)
    hrdThr = single(0.5);  % hard-threshold for shrinkage coefficient
    [U_rc, ~] = eig(ptsInRom*ptsInRom');
    ptsInRom = U_rc'*ptsInRom;  % transform coefficient
    shcInRom = 1-(stdGlb^2)*((ptsInRom.^(-1)).^(2));  % ^(-1) is faster than ^(-2)
    shcInRom = shcInRom.*(shcInRom>=hrdThr);
    ptsInRom = ptsInRom.*shcInRom;  % not necessary inverse transform
end


% match patches for reference patch
cntNbhInds = single(size(rptIndNbh, 2));
K = size(rptIndNbh, 1);
N = min(simNum, min(cntNbhInds));
simIndNbh = zeros(K, N, 'single');
for k = 1: K
    curInds = rptIndNbh(k, :);  % current inds in the neighborhood
    cdiInds = curInds(2:cntNbhInds);  % candidate inds of refInd
    refPat = ptsInRom(:, curInds(1));
    cdiPats = ptsInRom(:, cdiInds);
    
    dist = sum((cdiPats-repmat(refPat, 1, cntNbhInds-1)).^2, 1);  % cost most time when full color patSize
    [~, indc]= sort(dist);
    simIndNbh(k, :) = [curInds(1), cdiInds(indc(1:(N-1)))];
end
end


