function Jm = algStep3vc(Nm, stdGlb, patSize, nbhRad, simNum, balDiv, RmIndP, balChl, simIndsRec)
%ALGSTEP3VC Summary of this function goes here
%   Detailed explanation goes here


% initializing
chlDim = size(Nm, 4);
RomSize = [size(Nm, 1), size(Nm, 2), RmIndP(1, 2), chlDim];
RctSize = RomSize - [patSize(1), patSize(1), patSize(2), 1] + 1;
ptnInPag = RctSize(1)*RctSize(2);  % pats number in a page of Bm

precis = single(1e-12);
JmUpd = zeros(size(Nm), 'single');
JmPut = JmUpd;

isBnr = false;  % for MPCAv.m
Kp = size(RmIndP, 1);
for kp = 1: Kp
    Ridp = RmIndP(kp, :);
    Rom = Nm(:, :, Ridp(1):Ridp(2), :);
    
    ptsInRom = im2Pats(Rom, patSize, RctSize);
    
    [uptInRom, prsInRom] = MPCAvc(stdGlb, ptsInRom, simIndsRec(:, :, kp), patSize, simNum, balChl, balDiv, isBnr);
    
    [RomUpd, RomPut] = pats2Im(uptInRom, prsInRom, patSize, nbhRad, RomSize, RctSize, ptnInPag, precis);
    
    JmUpd(:, :, Ridp(1):Ridp(2), :) = JmUpd(:, :, Ridp(1):Ridp(2), :) + RomUpd;
    JmPut(:, :, Ridp(1):Ridp(2), :) = JmPut(:, :, Ridp(1):Ridp(2), :) + RomPut;
end


Jm = JmUpd./JmPut;
end