function [RomUpd, RomPut] = pats2Im(uptInRom, prsInRom, patSize, nbhRad, RomSize, RctSize, ptnInPag, precis)
%PATS2IM Summary of this function goes here
%   Weighted in page, average between pages.


% reshape uptInRom, prsInRom
uptInRom = reshape(permute(uptInRom, [2, 1, 3]), [], patSize(1)*patSize(1), patSize(2), RctSize(4));
prsInRom = reshape(permute(prsInRom, [2, 1, 3]), [], patSize(1)*patSize(1), patSize(2), RctSize(4));


% pats2Im
zroPag = zeros([RomSize(1), RomSize(2), patSize(2), RctSize(4)], 'single');
RomUpd = zeros(RomSize, 'single');
RomPut = RomUpd;  % precise used times
for p = 1: (1+2*nbhRad(2))
    uptInPag = uptInRom((1+ptnInPag*(p-1)):(ptnInPag*p), :, :, :);
    prsInPag = prsInRom((1+ptnInPag*(p-1)):(ptnInPag*p), :, :, :);
    
    % main part of pats2Im
    PagUpd = zroPag;
    PagPrs = zroPag;
    for c = 1: patSize(1)
        for r = 1: patSize(1)
            dim = r+(c-1)*patSize(1);
            rEnd = r+RctSize(1)-1;
            cEnd = c+RctSize(2)-1;
            PagUpd(r:rEnd, c:cEnd, :, :) = reshape(uptInPag(:, dim, :, :), RctSize(1), RctSize(2), patSize(2), RctSize(4)) + PagUpd(r:rEnd, c:cEnd, :, :);
            PagPrs(r:rEnd, c:cEnd, :, :) = reshape(prsInPag(:, dim, :, :), RctSize(1), RctSize(2), patSize(2), RctSize(4)) + PagPrs(r:rEnd, c:cEnd, :, :);
        end
    end
    
    
    % put into RomUpd
    PagPrs(PagPrs<precis) = precis;
    Page = PagUpd./PagPrs;
    pEnd = p+patSize(2)-1;
    RomUpd(:, :, p:pEnd, :) = RomUpd(:, :, p:pEnd, :) + Page;
    RomPut(:, :, p:pEnd, :) = RomPut(:, :, p:pEnd, :) + single(PagPrs>precis);  % average between pages
end
end


