function iscWithinFnirs(params, datafiles) 
% Calculate ISC within a group of subjects (each subject to average of
% others). Saves results as .mat file
%
% ARGUMENTS:
%  - params: struct specifying experiment params with fields -
%           subs: {1x18 cell}
%          crop: [tr1 tr2] % start timepoint and end timepoint
%      savedir: ''
%         name: 'CNvideo_CNgroup' %name of resulting nifti
%         type: 'within' %type of ISC (within/between)
%
% - datafiles: list of .mat datafiles to run isc over. Each .mat file
%      contains a struct with the following fields:
%           oxyData: [3000*47 double] %timepoints x channels
%


fprintf(['\n *** Calculating ISC (' params.type '-group): ' params.name '***\n']);

%% load all data
fprintf('calculating sum of all...\n')

for i = 1:length(datafiles)
    fprintf([num2str(i) '-']);
    
    % load subdata
    data = load(datafiles{i});
    
    % aggregate data
    if i == 1
        data_all = data.oxyData;
    else
        data_all(:,:,i) = data.oxyData;
    end
end
%

%% Calc sum of all
sum_all = sum(data_all,3);

sum_savename = fullfile(params.savedir, [params.name '_sumAll.mat']);
save(sum_savename, 'sum_all', 'params');

%% Do ISC
fprintf('\ndoing isc...\n');

corr_data = nan(size(data_all,2),length(datafiles));

for i = 1:length(datafiles)
    
    fprintf([num2str(i) '-']);
    
    % load subdata
    subj = data_all(:,:,i);
    
    % get avg others
    others = (sum_all - subj) / (length(params.subs)-1);
    
    % crop and zscore
    subj = zscore(subj(params.crop(1):params.crop(2),:),[],1);
    others = zscore(others(params.crop(1):params.crop(2),:), [],1);
    
    % correlate    
    
    for j = 1:size(subj,2)
            corr_data(j,i) = corr(subj(:,j),others(:,j)); % ISC-based
    end
    
end

fprintf('done! \n Making files ...');

isc = nanmean(corr_data,length(size(corr_data)));

%% Save
% Save ISC data
all_savename = fullfile(params.savedir, [params.name '_subISC.mat']);
save(all_savename, 'corr_data', 'params');

isc_savename = fullfile(params.savedir, [params.name '_ISC.mat']);
save(isc_savename, 'isc', 'params');
fprintf('done! \n');
