%% Analyze the stimulation with the photometry data
% 

basepath = pwd;
%photometry_file = dir(fullfile(basepath, '*N11_striatum.mat'));
%load(photometry_file.name);

sampling_rate = 130; % sampling rate of photometry set up - 130

color = 1; % 0/pink for striatum, 1/green for HPC

%% Plot photometry with behavior
%{
hold on
figure(1)
plot(photom_var.timestamps, photom_var.grabDA_z, 'g'); %plot(photometry(:,1), photometry(:,2), 'g');
ylabel('z score');
xlabel('time');
xlim([photom_var.timestamps(1), photom_var.timestamps(length(photom_var.timestamps(~isnan(photom_var.timestamps))))]);
for j = 1:length(rewarded_times)
    xline(rewarded_times(j, 1), '-b');
end
for k = 1:length(nonrewarded_times)
    xline(nonrewarded_times(k, 1), '-r');
end
hold off

%}


%% Average photometry around stim pulses

window = 10; % window of time around lick to average
samples = window*sampling_rate;

zscore_matrix = nan(length(photometryData.pulse_times), (samples*2)+1); 
photometryData.pulse_times = photometryData.pulse_times';
% convert to seconds!
photometryData.pulse_times = photometryData.pulse_times / 1000;

% average photometry data within a specified time window around stim
for j = 1:length(photometryData.pulse_times)
    curr_stim_time = photometryData.pulse_times(j);
    [~, stim_idx] = min(abs(photometryData.timestamps - photometryData.pulse_times(j)));
    start_idx = stim_idx - samples;
    end_idx = stim_idx + samples;
    if start_idx >= 1 && end_idx <= length(photometryData.timestamps)
        zscore_matrix(j, :) = photometryData.grabDA_z(start_idx:end_idx);
    end
end

%new_zscore_matrix = zscore_matrix;
%new_zscore_matrix(any(isnan(new_zscore_matrix), 2), :) = [];

zscore_matrix(any(isnan(zscore_matrix), 2), :) = [];
med_z_stim = median(zscore_matrix, 1); % median at each timepoint
time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(zscore_matrix);                          % Number of eExperimentsn In Data Set
avg_z_stim = mean(zscore_matrix, 1);              % Mean Of All Experiments At Each Value Of ,x 
stim_SEM = std(zscore_matrix, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
stim_CI95 = bsxfun(@times, stim_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw

% statistical significance of peak before lick
stim_baseline = zscore_matrix(:,1:1040);
avg_stim_baseline = mean(stim_baseline, 1);
mn = mean(avg_stim_baseline);
st_d = std(avg_stim_baseline);
sample_mn = min(med_z_stim(1170:1430));
deg_free = length(photometryData.pulse_times)-1;
t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));

%% Plot  average photometry level around stim pulses
if color == 0
    plot_color = [0.960784313725490, 0.152941176470588, 0.905882352941176]; % pink
    % striatum_pink = [0.960784313725490, 0.152941176470588, 0.690196078431373];
    % pink_conf = [0.968627450980392, 0.466666666666667, 0.701960784313725];
else
    plot_color = [0.031372549019608, 0.470588235294118, 0.149019607843137]; % green
end


figure('color','white');
ax = gca;
ax.FontSize = 15;
hold on
plot(time, med_z_stim, 'Color', plot_color, 'LineWidth', 2);
% plot(time, med_z_stim, 'g', 'LineWidth', 2);

fill([time,fliplr(time)], [(stim_CI95(1,:)+med_z_stim),fliplr((stim_CI95(2,:)+med_z_stim))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
% fill([time,fliplr(time)], [(stim_CI95(1,:)+med_z_stim),fliplr((stim_CI95(2,:)+med_z_stim))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around stims');
grid on;
hold off

% 'color', pink_conf

figure('color','white');
customGreen = [0.2039 0.7294 0.2039];
plot(photometryData.timestamps, photometryData.grabDA_z, 'Color', customGreen)
xlabel('time (s)');
ylabel('z-score');
title('GRAB-DA in hippocampus');








