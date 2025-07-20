

clc; clear; close all


outNum = 16;
P_pca = PROw(outNum);


pcf = P_pca{2, 3}(:, end);
pcf = sign(pcf(1))*pcf;
fprintf('%d outliers, first principal component [%.4f; %.4f]\n', outNum, pcf)








