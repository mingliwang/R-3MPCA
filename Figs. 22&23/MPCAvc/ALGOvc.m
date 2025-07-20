function [Jm, PSMA] = ALGOvc(Im, Nm, stdGau)
%ALGOVC Algorithm for color video
%   Detailed explanation goes here


[cimTfm, stdGlb] = clorTfm(stdGau, Nm);
Nm = permute(Nm, [1, 2, 4, 3]);
ImSize = size(Nm);
Nm = reshape(reshape(Nm, [], ImSize(4))*cimTfm, ImSize);


Pars = parsSet('video');
[rptNbhRC, RmIndP, rptNbhP] = refeNbh(ImSize(1:3), Pars.patSize, Pars.nbhRad, Pars.rptStep, Pars.cptStep);


Tm = algStep1vc(Nm, stdGlb, Pars.patSize, Pars.nbhRad, Pars.simNum, Pars.balDiv, ...
    RmIndP, rptNbhRC, rptNbhP);
Jm = permute(reshape(reshape(Tm, [], ImSize(4))/cimTfm, ImSize), [1, 2, 4, 3]);
PSMA(1, :) = [psnr(uint8(Jm), uint8(Im)), mean(double(abs(uint8(Jm(:))-uint8(Im(:)))))];


[Tm, simIndsRec] = algStep2vc(Nm, stdGlb, Pars.patSize, Pars.nbhRad, Pars.simNum, Pars.balDiv, ...
    RmIndP, rptNbhRC, rptNbhP, Tm, Pars.balCps, Pars.balChl);
Jm = permute(reshape(reshape(Tm, [], ImSize(4))/cimTfm, ImSize), [1, 2, 4, 3]);
PSMA(2, :) = [psnr(uint8(Jm), uint8(Im)), mean(double(abs(uint8(Jm(:))-uint8(Im(:)))))];


if (stdGlb>=15)
    Tm = algStep3vc(Tm, 0.1*stdGlb, Pars.patSize, Pars.nbhRad, Pars.simNum, Pars.balDiv, ...
        RmIndP, Pars.balChl, simIndsRec);
    Jm = permute(reshape(reshape(Tm, [], ImSize(4))/cimTfm, ImSize), [1, 2, 4, 3]);
    PSMA(3, :) = [psnr(uint8(Jm), uint8(Im)), mean(double(abs(uint8(Jm(:))-uint8(Im(:)))))];
end
end


