%% Analyze the behavioral data with the photometry data

% get behavior data
patch_behav = getPatchBehavior;

% load photometry data
data_folder = '/Users/lucyanderson/Library/Mobile Documents/com~apple~CloudDocs/Buzsaki/pyphotometry/N7/10.28/';
file = 'N7_striatum-2024-10-28-114829_photometry'; 
data_filename = strcat(file, '.mat');
data = strcat(data_folder, data_filename);
load(data);

rewarded_times = [];
nonrewarded_times = [];
% divide rewarded times vs nonrewarded times
for i = 1:length(patch_behav.timestamps)
    if patch_behav.reward_outcome(i) == 0
        nonrewarded_times = [nonrewarded_times; patch_behav.timestamps(i)];
    else
        rewarded_times = [rewarded_times; patch_behav.timestamps(i)];
    end
end


% synchronize
sync = getSyncPhotometry(photometryData);

hold on
figure
plot(sync(:,1), sync(:,2), 'g');
xlim([0, sync(length(sync),1)]);
for j = 1:length(rewarded_times)
    xline(rewarded_times(j), '-b');
end
for k = 1:length(nonrewarded_times)
    xline(nonrewarded_times(k), '-r');
end
hold off





