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



% pyphotometry sampling rate: 130
window = 1; % window of time around lick to average
sampling_rate = 130;
samples = window*sampling_rate;
zscore_matrix = nan(length(rewarded_times), (samples*2)+1);
zscore_matrix_non = nan(length(nonrewarded_times), (samples*2)+1);

for j = 1:length(rewarded_times)
    curr_reward_time = rewarded_times(j);
    [~, reward_idx] = min(abs(photometryData.timestamps - rewarded_times(j)));
   % photometryData.timestamps(sampling_rate*window)
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= size(photometryData.timestamps, 2)
        zscore_matrix(j, :) = sync(start_idx:end_idx, 2);

    end
end

avg_z_reward = mean(zscore_matrix, 1);
time = linspace(-window, window, ((samples*2)+1));

for k = 1:length(nonrewarded_times)
    curr_reward_time = nonrewarded_times(j);
    [~, reward_idx] = min(abs(photometryData.timestamps - nonrewarded_times(k)));
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= size(photometryData.timestamps, 2)
        zscore_matrix_non(k, :) = sync(start_idx:end_idx, 2);

    end
end

avg_z_no_reward = mean(zscore_matrix_non, 1);

% plot 
figure(2);
subplot(2,1,1)
plot(time, avg_z_reward, 'LineWidth', 2);
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Rewards');
grid on;

subplot(2,1,2)
plot(time, avg_z_no_reward, 'LineWidth', 2);
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Non-rewarded Licks');
grid on;

%% Average photometry across licks

window = 5; % window of time around lick to average
sampling_rate = 130; % sampling rate of photometry set up
samples = window*sampling_rate;

zscore_matrix = nan(length(rewarded_times), (samples*2)+1); % rewarded
zscore_matrix_non = nan(length(nonrewarded_times), (samples*2)+1); % not rewarded


% average photometry data within a specified time window around rewarded
% licks
for j = 1:length(rewarded_times)
    curr_reward_time = rewarded_times(j);
    [~, reward_idx] = min(abs(photometryData.timestamps - rewarded_times(j)));
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= size(photometryData.timestamps, 2)
        zscore_matrix(j, :) = sync(start_idx:end_idx, 2);

    end
end

med_z_reward = median(zscore_matrix, 1); % median at each timepoint
time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(zscore_matrix);                          % Number of eExperimentsn In Data Set
avg_z_reward = mean(zscore_matrix, 1);              % Mean Of All Experiments At Each Value Of ,x 
reward_SEM = std(zscore_matrix, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
reward_CI95 = bsxfun(@times, reward_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


% average photometry data within a specified time window around nonrewarded
% licks
for k = 1:length(nonrewarded_times)
    curr_reward_time = nonrewarded_times(j);
    [~, reward_idx] = min(abs(photometryData.timestamps - nonrewarded_times(k)));
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= size(photometryData.timestamps, 2)
        zscore_matrix_non(k, :) = sync(start_idx:end_idx, 2);

    end
end

med_z_no_reward = median(zscore_matrix_non, 1);

% calculate confidence intervals
N = height(zscore_matrix_non);                          % Number of experimentsn In Data Set
avg_z_no_reward = mean(zscore_matrix_non, 1);              % Mean Of All Experiments At Each Value Of ,x 
nonreward_SEM = std(zscore_matrix_non, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
nonreward_CI95 = bsxfun(@times, nonreward_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


%% Plot  average photometry level around licks
figure(2);
subplot(2,1,1)
plot(time, med_z_reward, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Rewards');
grid on;
hold off

subplot(2,1,2)
plot(time, med_z_no_reward, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(nonreward_CI95(1,:)+avg_z_no_reward),fliplr((nonreward_CI95(2,:)+avg_z_no_reward))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Non-rewarded Licks');
grid on;
hold off


