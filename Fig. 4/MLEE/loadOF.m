function Yr = loadOF(outNum)
%LOADOF Load the Old Faithful dataset
%   Detailed explanation goes here


load('dataset/OldFaithful.mat')
N = size(OldFaithful, 1);

% Normalized: mean=0, and let sep cov be 1
Y = double(OldFaithful');
Ym = mean(Y, 2);
Yr = Y - repmat(Ym, 1, N);
cov = (Yr*Yr')/N;
Yr(1, :) = Yr(1, :)/sqrt(cov(1, 1));
Yr(2, :) = Yr(2, :)/sqrt(cov(2, 2));
% plot(Yr(1, :), Yr(2, :), 'o')

% Outlier: add outlier
x_otl = 16;
outlier = [x_otl; -1];
randn('seed', 6);
nois = randn(size(Y, 1), outNum-1);
if (outNum==0)
    outliers = [];
else
    outliers = cat(2, outlier, repmat(outlier, 1, outNum-1) + 0.2*nois);
end

% Output data
Yr = double(cat(2, Yr, outliers));
end

