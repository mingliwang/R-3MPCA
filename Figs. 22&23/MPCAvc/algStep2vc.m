function [Jm, simIndNbhRec] = algStep2vc(Nm, stdGlb, patSize, nbhRad, simNum, balDiv, RmIndP, rptNbhRC, rptNbhP, Pm, balCps, balChl)
%ALGSTEP2VC Summary of this function goes here
%   Detailed explanation goes here


Pm = balCps*Pm+(1-balCps)*Nm;

% initializing
chlDim = size(Nm, 4);
RomSize = [size(Nm, 1), size(Nm, 2), RmIndP(1, 2), chlDim];
RctSize = RomSize - [patSize(1), patSize(1), patSize(2), 1] + 1;
ptnInPag = RctSize(1)*RctSize(2);  % pats number in a page of Bm

[Krc, Nrc] = size(rptNbhRC);
Np = size(rptNbhP, 2);
rptIndNbh0 = repmat(single(rptNbhRC), 1, 1, Np);
precis = single(1e-12);
JmUpd = zeros(size(Nm), 'single');
JmPut = JmUpd;

isBnr = true;
Kp = size(RmIndP, 1);
simIndNbhRec = zeros(Krc, simNum, Kp, 'single');
for kp = 1: Kp
    Ridp = RmIndP(kp, :);
    Rom = Nm(:, :, Ridp(1):Ridp(2), :);
    Pom = Pm(:, :, Ridp(1):Ridp(2), :);
    
    rptIndNbh = rptIndNbh0 + ptnInPag*repmat(reshape(rptNbhP(kp, :)-1, 1, 1, []), Krc, Nrc, 1);
    rptIndNbh = reshape(rptIndNbh, Krc, [], 1);  % first column is refPatInds
    
    ptsInRom = im2Pats(Rom, patSize, RctSize);
    
    simIndNbhRec(:, :, kp) = simiMat(0, simNum, im2Pats(Pom, patSize, RctSize), rptIndNbh);
    [uptInRom, prsInRom] = MPCAvc(stdGlb, ptsInRom, simIndNbhRec(:, :, kp), patSize, simNum, balChl, balDiv, isBnr);
    
    [RomUpd, RomPut] = pats2Im(uptInRom, prsInRom, patSize, nbhRad, RomSize, RctSize, ptnInPag, precis);
    
    JmUpd(:, :, Ridp(1):Ridp(2), :) = JmUpd(:, :, Ridp(1):Ridp(2), :) + RomUpd;
    JmPut(:, :, Ridp(1):Ridp(2), :) = JmPut(:, :, Ridp(1):Ridp(2), :) + RomPut;
end


Jm = JmUpd./JmPut;
end