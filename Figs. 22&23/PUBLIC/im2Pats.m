function ptsInRom = im2Pats(Rom, patSize, RctSize)
%IM2PATS  Summary of this function goes here
%   Merge pagsTrf.


% im2Pats
ptsInRom = zeros(prod(RctSize(1:3)), patSize(1)*patSize(1)*patSize(2), RctSize(4), 'single');
for p = 1: patSize(2)
    for c = 1: patSize(1)
        for r = 1: patSize(1)
            Rect = Rom(r:(r+RctSize(1)-1), c:(c+RctSize(2)-1), p:(p+RctSize(3)-1), :);
            coln = r+(c-1)*patSize(1)+(p-1)*(patSize(1)*patSize(1));
            ptsInRom(:, coln, :) = reshape(Rect, [], RctSize(4));
        end
    end
end
ptsInRom = permute(ptsInRom, [2, 1, 3]);  % for simiMat and MPCA
end


