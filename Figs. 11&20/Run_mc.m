

clc; clear; close all


figure('Position', [200, 200, 500, 300]);


X_orig = importdata('dataset/CHINADRESS_c3.mat');


Ours2 = PRMc_tvs(X_orig, 0);


