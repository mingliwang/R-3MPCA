function [Jm, PSMA] = ALGOv(Im, Nm, stdGau)
%ALGOV Algorithm for volumetric data
%   Detailed explanation goes here


Pars = parsSet('vlum', stdGau);

[rptNbhRC, RmIndP, rptNbhP] = refeNbh(cat(2, size(Nm, 1), size(Nm, 2), size(Nm, 3)), ...
    Pars.patSize, Pars.nbhRad, Pars.rptStep, Pars.cptStep);

Jm = algStep1v(Nm, stdGau, Pars.patSize, Pars.nbhRad, Pars.simNum, Pars.balDiv, ...
    RmIndP, rptNbhRC, rptNbhP, true);
PSMA(1, :) = [psnr(uint8(Jm), uint8(Im)), mean(double(abs(uint8(Jm(:))-uint8(Im(:)))))];

if (stdGau>=15)
    [Jm, simIndRec] = algStep2v(Nm, stdGau, Pars.patSize, Pars.nbhRad, Pars.simNum, ...
        RmIndP, rptNbhRC, rptNbhP, Jm, Pars.balCps, Pars.balChl, Pars.balWie);
    PSMA(2, :) = [psnr(uint8(Jm), uint8(Im)), mean(double(abs(uint8(Jm(:))-uint8(Im(:)))))];
    
    Jm = algStep3v(Jm, 0.1*stdGau, Pars.patSize, Pars.nbhRad, Pars.simNum, Pars.balDiv, ...
        RmIndP, false, simIndRec);
    PSMA(3, :) = [psnr(uint8(Jm), uint8(Im)), mean(double(abs(uint8(Jm(:))-uint8(Im(:)))))];
elseif (stdGau>=10)
    Jm = algStep1v(Jm, 0.1*stdGau, Pars.patSize, Pars.nbhRad, Pars.simNum, Pars.balDiv, ...
        RmIndP, rptNbhRC, rptNbhP, false);
    PSMA(2, :) = [psnr(uint8(Jm), uint8(Im)), mean(double(abs(uint8(Jm(:))-uint8(Im(:)))))];
end
end


