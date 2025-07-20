

clc; clear; close all


outNum = 81;
Pw_pca = PROw_L1(outNum);

pcf = Pw_pca{2, 3}(:, end);
pcf = sign(pcf(1))*pcf;

fprintf('warm start, %d outliers, first principal component [%.4f; %.4f]\n', outNum, pcf)








