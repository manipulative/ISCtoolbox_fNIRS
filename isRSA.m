function [rho p] = isRSA(data1, data2)
% function used to calculate similartiy matrix based on matrix data
% Input:
%      data1 & data2: two similarity matrix, [subject x subject]
% Output
%      [rho p]: spearman correlation and p value between two matrix

[a,b] = size(data1);
d1_vec = data1(tril(ones(a,b),-1)>0);
[a,b] = size(data2);
d2_vec = data2(tril(ones(a,b),-1)>0);

[rho p] = corr(d1_vec, d2_vec, 'Type','Spearman');

