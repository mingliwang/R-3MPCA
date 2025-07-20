function [rptNbhRC, RmIndP, rptNbhP] = refeNbh(ImSize, patSize, nbhRad, rptStep, cptStep)
%REFENBH Summary of this function goes here
%   Detailed explanation goes here


rptNbhRC = refeNbhRC(ImSize(1:2), patSize(1), nbhRad(1), rptStep(1), cptStep(1));

[RmIndP, rptNbhP] = refeNbhP(ImSize(3), patSize(2), nbhRad(2), rptStep(2), cptStep(2));
end


function rptNbh = refeNbhRC(ImSize, patSize, nbhRad, rptStep, cptStep)
%REFENBHRC refPat's cdiPats inds, row-column direction
%   Detailed explanation goes here


% refPat in (row, column) direction, search window inds
[rptIndSwR, rswIndR, nbhRadR] = refeInds(ImSize(1), patSize, nbhRad, rptStep);
[rptIndSwC, rswIndC, nbhRadC] = refeInds(ImSize(2), patSize, nbhRad, rptStep);

% cdiPat with cptStep, in search window
cptIndSwR = single(0:cptStep:nbhRadR)';
if (cptIndSwR(end)~=nbhRadR)
    cptIndSwR(end+1) = nbhRadR;
end
cptIndSwR = cat(1, flip(-cptIndSwR(2:end), 1), cptIndSwR) + (nbhRadR+1);

cptIndSwC = single(0:cptStep:nbhRadC)';
if (cptIndSwC(end)~=nbhRadC)
    cptIndSwC(end+1) = nbhRadC;
end
cptIndSwC = cat(1, flip(-cptIndSwC(2:end), 1), cptIndSwC) + (nbhRadC+1);

cptNumSw = [length(cptIndSwR), length(cptIndSwC)];
cptIndSw = zeros(prod(cptNumSw), 1, 'single');
cnt = 1;
for j = 1: cptNumSw(2)
    for i = 1: cptNumSw(1)
        cptIndSw(cnt) = cptIndSwR(i) + (cptIndSwC(j)-1)*(2*nbhRadR+1);
        cnt = cnt + 1;
    end
end

% refPat's cdiPats inds, row-column direction
rectSize = [ImSize(1)-patSize+1, ImSize(2)-patSize+1];
patIndMap = reshape(single(1:1:(rectSize(1)*rectSize(2)))', rectSize(1), rectSize(2));

rptNum = [size(rptIndSwR, 1), size(rptIndSwC, 1)];
rptInd = zeros(prod(rptNum), 1, 'single');
rptNbh = zeros((2*nbhRadR+1)*(2*nbhRadC+1), prod(rptNum), 'single');
isCptStpOne = (cptStep==1);
cdiIndSw0 = single(1:1:((2*nbhRadR+1)*(2*nbhRadC+1)));  % for the case ~isCptStpOne
cnt = 1;
for j = 1: rptNum(2)
    for i = 1: rptNum(1)
        rptIndSw = rptIndSwR(i) + (rptIndSwC(j)-1)*(2*nbhRadR+1);
        tmp = patIndMap(rswIndR(i, 1):(rswIndR(i, 2)-(patSize-1)), rswIndC(j, 1):(rswIndC(j, 2)-(patSize-1)));
        tmp = tmp(:);
        rptInd(cnt) = tmp(rptIndSw);
        if isCptStpOne
            tmp(rptIndSw) = [];
            rptNbh(:, cnt) = [rptInd(cnt); tmp];
        else
            cdiIndSw = cdiIndSw0;
            cdiIndSw([rptIndSw; cptIndSw]) = [];
            tmp(rptIndSw) = tmp(cdiIndSw(1));  % replace tmp(rptIndSw) with the other one
            rptNbh(:, cnt) = tmp;
        end
        cnt = cnt + 1;
    end
end
if isCptStpOne
    rptNbh = rptNbh';
else
    rptNbh = cat(2, rptInd, rptNbh(cptIndSw, :)');
end
end


function [RmIndP, rptNbhRm] = refeNbhP(ImSizeP, patSize, nbhRad, rptStep, cptStep)
%REFENBHP refPat's inds, page direction
%   Detailed explanation goes here


% refPat inds in Rom, refRom inds in Rom, page direction
[rptIndRm, RmIndP, nbhRadP] = refeInds(ImSizeP, patSize, nbhRad, rptStep);

% cdiPat with cptStep, in refRom, page direction
cptIndRm = single(0:cptStep:nbhRadP);
if (cptIndRm(end)~=nbhRadP)
    cptIndRm(end+1) = nbhRadP;
end
cptIndRm = cat(2, flip(-cptIndRm(2:end), 2), cptIndRm) + (nbhRadP+1);

% refPat's cdiPats inds in Rm (locally), page direction
rptNumP = length(rptIndRm);
rptNbhRm = zeros(rptNumP, 2*nbhRadP+1, 'single');

isCptStpOne = (cptStep==1);
patIndRm = single(1:1:(2*nbhRadP+1));
for k = 1: rptNumP
    tmp = patIndRm;
    if isCptStpOne
        tmp(rptIndRm(k)) = [];
        rptNbhRm(k, :) = [rptIndRm(k), tmp];
    else
        cdiIndRm = patIndRm;
        cdiIndRm([rptIndRm(k), cptIndRm]) = [];
        tmp(rptIndRm(k)) = tmp(cdiIndRm(1));  % replace tmp(rptIndSw) with the other one
        rptNbhRm(k, :) = tmp;
    end
end
if ~isCptStpOne
    rptNbhRm = cat(2, rptIndRm, rptNbhRm(:, cptIndRm));
end
end


function [rptInSwSub, rswIndSub, nbhRadSub] = refeInds(ImSizeSub, patSizeSub, nbhRadSub, rptStepSub)
%REFEINDS Summary of this function goes here
%   Detailed explanation goes here


% refPat inds, in a (row, column, page) direction
rptIndSub = single((1:rptStepSub:(ImSizeSub-patSizeSub+1)))';
if (rptIndSub(end)~=ImSizeSub-patSizeSub+1)
    rptIndSub(end+1) = ImSizeSub-patSizeSub+1;  % to the end
end
rptIndSubNum = length(rptIndSub);

% refPat inds in search window, refSw, in a (row, column, page) direction
nbhRadSub = min(nbhRadSub, floor((ImSizeSub-patSizeSub)/2));  % search window radius
rswIndSub = [rptIndSub-nbhRadSub, rptIndSub+(patSizeSub-1)+nbhRadSub];
outsideSub = [1-rswIndSub(:, 1), rswIndSub(:, 2)-ImSizeSub];

rptInSwSub = (nbhRadSub+1)*ones(rptIndSubNum, 1, 'single');
rptInSwSub = rptInSwSub - outsideSub(:, 1).*single(outsideSub(:, 1)>0);
rptInSwSub = rptInSwSub + outsideSub(:, 2).*single(outsideSub(:, 2)>0);  % rpt inds in search window

rswIndSub = rswIndSub + repmat(outsideSub(:, 1).*single(outsideSub(:, 1)>0), 1, 2);
rswIndSub = rswIndSub - repmat(outsideSub(:, 2).*single(outsideSub(:, 2)>0), 1, 2);
end


