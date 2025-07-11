%% Analyze the sleep photometry data 
% 
% Takes in photometry data from the pre- and post- sleep sessions that has 
% not already been synced to intan. Plots the averaged photometry data around 
% ripple events for the entire session.

color = 0; % 0/pink for striatum, 1/green for HPC

%% combine all photometry sessions
% full_photometry.timestamps = sleep_sync.timestamps;
% full_photometry.grabDA_z = sleep_sync.grabDA_z;
% 
% full_photometry.timestamps = [full_photometry.timestamps; photometry.timestamps];
% full_photometry.grabDA_z = [full_photometry.grabDA_z; photometry.grabDA_z];
% 
% full_photometry.timestamps = [full_photometry.timestamps; sleep_sync.timestamps+MergePoints.timestamps(3,1)];
% full_photometry.grabDA_z = [full_photometry.grabDA_z; sleep_sync.grabDA_z];
% 
% figure
% plot(full_photometry.timestamps, full_photometry.grabDA_z);


%% Load variables

% set pre_post to be 1 for pre behavior sleep, 2 for post behavior sleep  
pre_post = 1;

basepath = pwd;
[~, currentFolderName] = fileparts(basepath);

sampling_rate = 130; % sampling rate of photometry set up - 130

if isempty(dir(fullfile(basepath, '*SleepPhotomSynced.mat')))
    photometry_file = dir(fullfile(basepath, '*_photometry.mat'));
    load(photometry_file.name);

    % synchronize time stamps for sleep session
    sleep_sync = getSyncPhotometry(photometryData);

    C = strsplit(pwd,'\');
    save([basepath filesep C{end} '.SleepPhotomSynced.mat'],'sleep_sync');
else
    photometry_file = dir(fullfile(basepath, '*SleepPhotomSynced.mat'));
    load(photometry_file.name);
end

ripple_file = dir(fullfile(basepath, '..', '*ripples.events.mat'));
load(ripple_file.name)

merge_pt_file = dir(fullfile(basepath, '..', '*MergePoints.events.mat'));
load(merge_pt_file.name)

if pre_post == 1 
    sleep_start = MergePoints.timestamps(1,1); % start of sleep session
    sleep_end = MergePoints.timestamps(1,2); % time stamp of end of sleep session 
elseif pre_post == 2
    sleep_start = MergePoints.timestamps(3,1); 
    sleep_end = MergePoints.timestamps(3,2);  
end

% make this a fxn
%function [avg] = avgAcrossEvents(event)

%{
    sleep_start = MergePoints.timestamps(2,1); 
    sleep_end = MergePoints.timestamps(2,2); 
%}

color = 0; % 0/pink for striatum, 1/green for HPC

%% combine all photometry sessions
% full_photometry.timestamps = sleep_sync.timestamps;
% full_photometry.grabDA_z = sleep_sync.grabDA_z;
% 
% full_photometry.timestamps = [full_photometry.timestamps; photometry.timestamps];
% full_photometry.grabDA_z = [full_photometry.grabDA_z; photometry.grabDA_z];
% 
% full_photometry.timestamps = [full_photometry.timestamps; sleep_sync.timestamps+MergePoints.timestamps(3,1)];
% full_photometry.grabDA_z = [full_photometry.grabDA_z; sleep_sync.grabDA_z];


%% full trace

% BLUE = NREM
% RED = REM
% YELLOW = WAKE

if exist('full_photometry')
    sleep_state_file = dir(fullfile(basepath, '..', '*SleepState.states.mat'));
    load(sleep_state_file.name)
    adjusted_ts = sleep_sync.timestamps + sleep_start;
    
    darkerGreen = [0.1098    0.6000    0.2392];
    striatum_pink = [0.960784313725490, 0.152941176470588, 0.905882352941176];
    figure('color','white')
    ax = gca;
    ax.FontSize = 15;
    hold on
    plot(adjusted_ts, sleep_sync.grabDA_z, 'Color', darkerGreen, 'LineWidth', 2);
    for ii = 1:size(SleepState.ints.NREMstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.NREMstate(ii, 1), SleepState.ints.NREMstate(ii, 2), ...
            SleepState.ints.NREMstate(ii, 2), SleepState.ints.NREMstate(ii, 1)];
        patch(x, y, [0.301960784313725   0.745098039215686   0.933333333333333], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    
    for ii = 1:size(SleepState.ints.REMstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.REMstate(ii, 1), SleepState.ints.REMstate(ii, 2), ...
            SleepState.ints.REMstate(ii, 2), SleepState.ints.REMstate(ii, 1)];
        patch(x, y, [1.000000000000000   0.266666666666667                   0], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    
    for ii = 1:size(SleepState.ints.WAKEstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.WAKEstate(ii, 1), SleepState.ints.WAKEstate(ii, 2), ...
            SleepState.ints.WAKEstate(ii, 2), SleepState.ints.WAKEstate(ii, 1)];
        patch(x, y, [1.000000000000000   1.000000000000000   0.066666666666667], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    xlim([adjusted_ts(1), adjusted_ts(end)]);
    ylim([min(y), max(y)])
    hold off

    sleep_state_file = dir(fullfile(basepath, '..', '*SleepState.states.mat'));
    load(sleep_state_file.name)
    adjusted_ts = sleep_sync.timestamps + sleep_start;
    
    darkerGreen = [0.1098    0.6000    0.2392];
    striatum_pink = [0.960784313725490, 0.152941176470588, 0.905882352941176];
    figure('color','white')
    ax = gca;
    ax.FontSize = 15;
    hold on
    plot(adjusted_ts, sleep_sync.grabDA_z, 'Color', darkerGreen, 'LineWidth', 2);
    for ii = 1:size(SleepState.ints.NREMstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.NREMstate(ii, 1), SleepState.ints.NREMstate(ii, 2), ...
            SleepState.ints.NREMstate(ii, 2), SleepState.ints.NREMstate(ii, 1)];
        patch(x, y, [0.301960784313725   0.745098039215686   0.933333333333333], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    
    for ii = 1:size(SleepState.ints.REMstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.REMstate(ii, 1), SleepState.ints.REMstate(ii, 2), ...
            SleepState.ints.REMstate(ii, 2), SleepState.ints.REMstate(ii, 1)];
        patch(x, y, [1.000000000000000   0.266666666666667                   0], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    
    for ii = 1:size(SleepState.ints.WAKEstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.WAKEstate(ii, 1), SleepState.ints.WAKEstate(ii, 2), ...
            SleepState.ints.WAKEstate(ii, 2), SleepState.ints.WAKEstate(ii, 1)];
        patch(x, y, [1.000000000000000   1.000000000000000   0.066666666666667], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    xlim([adjusted_ts(1), adjusted_ts(end)]);
    ylim([min(y), max(y)])
    hold off
end


%% Analyze sleep photometry

ripple_period = ripples.peaks(ripples.peaks <= sleep_end & ripples.peaks >= sleep_start);

% find average time between ripple events
ripple_durs = zeros(length(ripple_period), 1);
for k = 2:length(ripple_period)
    ripple_durs(k) = ripple_period(k)-ripple_period(k-1);
end
median(ripple_durs)

window = 5; % window of time around ripple to average
samples = window*sampling_rate;

% initialize matrix
ripple_matrix = nan(length(ripple_period), (samples*2)+1); 


% average photometry data within a specified time window around ripples
adjusted_ts = sleep_sync.timestamps + sleep_start;
for j = 1:length(ripple_period)
    [~, ripple_idx] = min(abs(adjusted_ts - ripple_period(j)));
    start_idx = ripple_idx - samples;
    end_idx = ripple_idx + samples;
    if start_idx >= 1 && end_idx <= height(adjusted_ts)
        ripple_matrix(j, :) = sleep_sync.grabDA_z(start_idx:end_idx);
        
        % % z score
        % baseline = ripple_matrix(j, :);
        % 
        % % Calculate baseline mean and std
        % baseline_mean = mean(baseline);
        % baseline_std = std(baseline);
        % 
        % % Z-score the entire trial using the baseline stats
        % ripple_matrix(j, :) = (ripple_matrix(j, :) - baseline_mean) / baseline_std;
    end
end

ripple_matrix(any(isnan(ripple_matrix), 2), :) = [];

% ripple_matrix_2 = ripple_matrix;
% new = [ripple_matrix_2; ripple_matrix];
% ripple_matrix = new;

median_ripple = median(ripple_matrix, 1); % median at each timepoint
time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(ripple_matrix);                          % Number of Experiments In Data Set
avg_ripple = mean(ripple_matrix, 1);              % Mean Of All Experiments At Each Value Of ,x 
ripple_SEM = std(ripple_matrix, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw
smooth_CI95 = smoothdata(ripple_CI95, 2);

% statistical significance of peak before lick
ripple_baseline = ripple_matrix(:,[1:120]);
avg_ripple_baseline = mean(ripple_baseline, 1);
mn = mean(avg_ripple_baseline);
st_d = std(avg_ripple_baseline);
sample_mn = min(median_ripple);
deg_free = length(ripple_period)-1;
t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));

%% Plot photometry against ripples

if color == 0
    % avg_color = [ 0.2392    0.2863    0.9608];
    % conf_color = 'b';
    avg_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
    conf_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
else
    avg_color = 'g';
    conf_color = [0.7176    0.9412    0.1020];
end

figure('color','white');
plot(time, median_ripple, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Ripples');% (Pre-task Sleep)');
grid on;
hold off


smoothed = smoothdata(median_ripple);
figure('color','white');
plot(time, smoothed, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(smooth_CI95(1,:)+smoothed),fliplr((smooth_CI95(2,:)+smoothed))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Ripples (smoothed)');
grid on;
hold off

%{
figure('color','white');
plot(time, median_ripple, 'color', avg_color, 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Average Z-score Around Ripples - ', currentFolderName, ' HPC'], [num2str(size(ripple_matrix, 1)), ' ripple events']);% (Pre-task Sleep)');
grid on;
hold off
%}

% chan = sessionInfo.AnatGrps(5).Channels;
% lfp = bz_GetLFP(chan);
% bz_eventCSD(lfp, ripples.peaks);

%{

%photom_pre_sleep = photometryData;
% photometry across sleep and behav
figure('color','white');
subplot(2,1,1)
plot(photom_pre_sleep.timestamps, photom_pre_sleep.grabDA_z)
%xlim([0, photom_pre_sleep.timestamps(size(photom_pre_sleep.timestamps))]);
ylabel('z score');
title('sleep DA');

subplot(2,1,2)
plot(photometryData.timestamps, photometryData.grabDA_z)
%xlim([0, photometryData.timestamps(size(photometryData.timestamps))]);
ylabel('z score');
title('behavior DA');


%}




