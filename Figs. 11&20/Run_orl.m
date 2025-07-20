

clc; clear; close all


figure('Position', [200, 200, 500, 300]);

X_orig = importdata('dataset/ORL_35.mat');

ind = 5;

Ours3 = PROrl_tw(X_orig(:, :, ind), ind-1);

Ours2b = PROrl_tvs(X_orig, 0);

Ours2a = PROrl_tvs(X_orig(:, :, ind), ind-1);


