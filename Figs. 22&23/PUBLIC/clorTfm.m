function [cimTfm, stdGlb] = clorTfm(stdGau, Jm)
%CLORTFM Color image (space) transformation
%

% noise variance stabilizatioin
stdGlb = min(single(stdGau));  % define global standard deviation
if (sum(stdGau)==0)
    stdTfm = single(eye(3));
else
    stdTfm = stdGlb*diag(stdGau)^(-1);
end

% color space transformation matrix; opponent color space or PCA
if max(stdGau)-stdGlb <= 2
    clrTfm = single([1/3 1/3 1/3; 0.5 0 -0.5; 0.25 -0.5 0.25]);  % opponent transformation
    clrTfm = (clrTfm./repmat(sqrt(sum(clrTfm.^2, 2)), 1, size(clrTfm, 2)))';
else
    [~, ~, D, ~] = size(Jm);
    Jc = reshape(permute(Jm, [1, 2, 4, 3]), [], D)*stdTfm;  % variance stabilization
    [clrTfm, ~] = eig(Jc'*Jc);  % the correlation matrix, PCA
    clrTfm = flip(clrTfm, 2);
end

% color image transformation matrix
cimTfm = single(stdTfm*clrTfm);
end
