function Jm = algStep1v(Nm, stdGlb, patSize, nbhRad, simNum, balDiv, RmIndP, rptNbhRC, rptNbhP, isBnr)
%ALGSTEP1V Summary of this function goes here
%   Detailed explanation goes here


% initializing
RomSize = [size(Nm, 1), size(Nm, 2), RmIndP(1, 2), 1];
RctSize = RomSize - [patSize(1), patSize(1), patSize(2), 1] + 1;
ptnInPag = RctSize(1)*RctSize(2);  % pats number in a page of Rom

[Krc, Nrc] = size(rptNbhRC);
Np = size(rptNbhP, 2);
rptIndNbh0 = repmat(single(rptNbhRC), 1, 1, Np);
precis = single(1e-12);
JmUpd = zeros(size(Nm), 'single');
JmPut = JmUpd;

Kp = size(RmIndP, 1);
for kp = 1: Kp
    Ridp = RmIndP(kp, :);
    Rom = Nm(:, :, Ridp(1):Ridp(2));
    
    rptIndNbh = rptIndNbh0 + ptnInPag*repmat(reshape(rptNbhP(kp, :)-1, 1, 1, []), Krc, Nrc, 1);
    rptIndNbh = reshape(rptIndNbh, Krc, [], 1);  % first column is refPatInds
    
    ptsInRom = im2Pats(Rom, patSize, RctSize);
    
    simIndNbh = simiMat(stdGlb, simNum, ptsInRom, rptIndNbh);
    [uptInRom, prsInRom] = MPCAv(stdGlb, ptsInRom, simIndNbh, patSize, simNum, balDiv, isBnr);
    
    [RomUpd, RomPut] = pats2Im(uptInRom, prsInRom, patSize, nbhRad, RomSize, RctSize, ptnInPag, precis);
    
    JmUpd(:, :, Ridp(1):Ridp(2)) = JmUpd(:, :, Ridp(1):Ridp(2)) + RomUpd;
    JmPut(:, :, Ridp(1):Ridp(2)) = JmPut(:, :, Ridp(1):Ridp(2)) + RomPut;
end


Jm = JmUpd./JmPut;
end