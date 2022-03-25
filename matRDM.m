function r = matRDM(data, method)
% function used to calculate similartiy matrix based on matrix data
% Input:
%      data: matrix, [subject x item/timepoint]
%      method: 'pcorr': Pearson correlation
%              'spearman': Spearman correlation
%              'mdist': difference between mean of whole varible
%              'amdist': absolute difference between mean of whole varible
%              'euclidean': Euclidean distance at item-wise
%              'seuclidean': standrized Euclidean distance at item-wose
%              'cosine': cosine distance at item-wise
% Output
%      r: similarity matrix

T = data;
r = NaN(size(T,1),1);

for ii = 1:size(T,1)
    for jj = 1:size(T,1)
        
        subj1 = T(ii,:);
        subj2 = T(jj,:);
        
        if strcmp(method,'pcorr')
            r(ii,jj) = corr(subj1',subj2');
        elseif strcmp(method,'spearman')
            r(ii,jj) = corr(subj1',subj2', 'Type', 'Spearman');
        elseif strcmp(method,'mdist')
            r(ii,jj) = mean(subj1) - mean(subj2);
        elseif strcmp(method,'amdist')
            r(ii,jj) = abs(mean(subj1) - mean(subj2));
        else
            r(ii,jj) = pdist2(subj1,subj2,method);
        end
    end
end