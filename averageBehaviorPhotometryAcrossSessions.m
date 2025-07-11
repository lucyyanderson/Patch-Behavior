<<<<<<< HEAD
=======
% z score across behavior sessions
%
% USAGE
%    z score photometry data from across behavior sessions and plot the 
%    dopamine traces around rewarded or nonrewarded licks, in the HPC and striatum
% 
%
% INPUTS 
%    need to load in the z score matrices around rewarded licks and
%    nonrewarded licks, and then update code to reflect proper session
%    numbers
%
%    =========================================================================

%{
%% combine sessions

window = 5; % seconds
sampling_rate = 130;
samples = window * sampling_rate;
baseline_idx = 1:325;

time = linspace(-window, window, ((samples*2)+1));

% Preallocate matrices for normalized data
norm_session1 = zeros(size(zscore_matrix_18));
norm_session2 = zeros(size(zscore_matrix_19));
norm_session3 = zeros(size(zscore_matrix_22));

% Normalize each trial in each session
sessions = {zscore_matrix_18, zscore_matrix_19, zscore_matrix_22};
norm_sessions_reward = {norm_session1, norm_session2, norm_session3};

for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, baseline_idx);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions_reward{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end

% Extract the normalized sessions
% norm_session1 = norm_sessions{1};
% norm_session2 = norm_sessions{2};
% norm_session3 = norm_sessions{3};


%% non rewarded
% Preallocate matrices for normalized data
norm_session1 = zeros(size(zscore_matrix_non_18));
norm_session2 = zeros(size(zscore_matrix_non_19));
norm_session3 = zeros(size(zscore_matrix_non_22));

% Normalize each trial in each session
sessions = {zscore_matrix_non_18, zscore_matrix_non_19, zscore_matrix_non_22};
norm_sessions_nonreward = {norm_session1, norm_session2, norm_session3};

for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, baseline_idx);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions_nonreward{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end



%% PLOT

figure('color','white');
subplot(2, 1, 1)
hold on
ax = gca;
ax.FontSize = 15;

%green_map = 
% rewarded
for ii=1:length(norm_sessions_reward)
    current_norm = norm_sessions_reward{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = summer(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Average Z-score Around Rewards');
end
grid on;
hold off

% nonrewarded
subplot(2, 1, 2)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_nonreward)
    current_norm = norm_sessions_nonreward{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = summer(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Average Z-score Around non rewarded licks');
end
grid on;
hold off



%% striatum

baseline_idx = 1:520;

% Preallocate matrices for normalized data
norm_session1 = zeros(size(zscore_matrix_16));
norm_session2 = zeros(size(zscore_matrix_17));
norm_session3 = zeros(size(zscore_matrix_20));

% Normalize each trial in each session
sessions = {zscore_matrix_16, zscore_matrix_17, zscore_matrix_20};
norm_sessions_reward = {norm_session1, norm_session2, norm_session3};

for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, baseline_idx);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions_reward{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end



%% non rewarded
% Preallocate matrices for normalized data
norm_session1 = zeros(size(zscore_matrix_non_16));
norm_session2 = zeros(size(zscore_matrix_non_17));
norm_session3 = zeros(size(zscore_matrix_non_20));

% Normalize each trial in each session
sessions = {zscore_matrix_non_16, zscore_matrix_non_17, zscore_matrix_non_20};
norm_sessions_nonreward = {norm_session1, norm_session2, norm_session3};

for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, baseline_idx);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions_nonreward{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end



%% PLOT

figure('color','white');
subplot(2, 1, 1)
hold on
ax = gca;
ax.FontSize = 15;

% rewarded
for ii=1:length(norm_sessions_reward)
    current_norm = norm_sessions_reward{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = spring(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Average Z-score Around Rewards');
end
grid on;
hold off

% nonrewarded
subplot(2, 1, 2)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_nonreward)
    current_norm = norm_sessions_nonreward{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = spring(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Average Z-score Around non rewarded licks');
end
grid on;
hold off


%{
global_mean = mean(all_reward(:), baseline));
global_std = std(all_reward(:, baseline));

norm_reward = (all_reward-global_mean)/global_std;

med_z_reward = median(norm_reward(1:119, 1)); % median at each timepoint


all_reward=[zscore_matrix_18; zscore_matrix_19; zscore_matrix_22];
all_nonreward=[zscore_matrix_non_18; zscore_matrix_non_19; zscore_matrix_non_22];

%}

%}























>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

% z score across behavior sessions
%
% USAGE
%    z score photometry data from across behavior sessions and plot the 
%    dopamine traces around rewarded or nonrewarded licks, in the HPC and striatum
% 
%
% INPUTS 
%    need to load in the z score matrices around rewarded licks and
%    nonrewarded licks, and then update code to reflect proper session
%    numbers
%
%    =========================================================================


<<<<<<< HEAD
%% combine sessions - rewarded
=======
%% combine sessions
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

window = 5; % seconds
sampling_rate = 130;
samples = window * sampling_rate;
baseline_idx = 1:325;

time = linspace(-window, window, ((samples*2)+1));

% Preallocate matrices for normalized data
norm_session1 = zeros(size(zscore_matrix_18));
norm_session2 = zeros(size(zscore_matrix_19));
norm_session3 = zeros(size(zscore_matrix_22));
<<<<<<< HEAD
norm_session4 = zeros(size(zscore_matrix_24));


% Normalize each trial in each session
sessions = {zscore_matrix_18, zscore_matrix_19, zscore_matrix_22, zscore_matrix_24};
norm_sessions_reward = {norm_session1, norm_session2, norm_session3, norm_session4};
=======

% Normalize each trial in each session
sessions = {zscore_matrix_18, zscore_matrix_19, zscore_matrix_22};
norm_sessions_reward = {norm_session1, norm_session2, norm_session3};
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, :);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions_reward{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end

% Extract the normalized sessions
% norm_session1 = norm_sessions{1};
% norm_session2 = norm_sessions{2};
% norm_session3 = norm_sessions{3};


%% non rewarded
% Preallocate matrices for normalized data
norm_session1 = zeros(size(zscore_matrix_non_18));
norm_session2 = zeros(size(zscore_matrix_non_19));
norm_session3 = zeros(size(zscore_matrix_non_22));
<<<<<<< HEAD
norm_session4 = zeros(size(zscore_matrix_non_24));

% Normalize each trial in each session
sessions = {zscore_matrix_non_18, zscore_matrix_non_19, zscore_matrix_non_22, zscore_matrix_non_24};
norm_sessions_nonreward = {norm_session1, norm_session2, norm_session3, norm_session4};
=======

% Normalize each trial in each session
sessions = {zscore_matrix_non_18, zscore_matrix_non_19, zscore_matrix_non_22};
norm_sessions_nonreward = {norm_session1, norm_session2, norm_session3};
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, :);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions_nonreward{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end



<<<<<<< HEAD
%% PLOT HPC rewarded vs nonrewarded
=======
%% PLOT
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

figure('color','white');
subplot(2, 1, 1)
hold on
ax = gca;
ax.FontSize = 15;

%green_map = 
% rewarded
for ii=1:length(norm_sessions_reward)
    current_norm = norm_sessions_reward{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

<<<<<<< HEAD
    cmap = summer(5);
=======
    cmap = summer(4);
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Average Z-score Around Rewards');
end
grid on;
hold off

% nonrewarded
subplot(2, 1, 2)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_nonreward)
    current_norm = norm_sessions_nonreward{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

<<<<<<< HEAD
    cmap = summer(5);
=======
    cmap = summer(4);
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Average Z-score Around non rewarded licks');
end
grid on;
hold off



%% striatum

baseline_idx = 1:520;

% Preallocate matrices for normalized data
norm_session1 = zeros(size(zscore_matrix_16));
norm_session2 = zeros(size(zscore_matrix_17));
norm_session3 = zeros(size(zscore_matrix_20));
<<<<<<< HEAD
norm_session4 = zeros(size(zscore_matrix_25));

% Normalize each trial in each session
sessions = {zscore_matrix_16, zscore_matrix_17, zscore_matrix_20, zscore_matrix_25};
norm_sessions_reward = {norm_session1, norm_session2, norm_session3, norm_session4};
=======

% Normalize each trial in each session
sessions = {zscore_matrix_16, zscore_matrix_17, zscore_matrix_20};
norm_sessions_reward = {norm_session1, norm_session2, norm_session3};
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, :);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions_reward{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end



%% non rewarded
% Preallocate matrices for normalized data
norm_session1 = zeros(size(zscore_matrix_non_16));
norm_session2 = zeros(size(zscore_matrix_non_17));
norm_session3 = zeros(size(zscore_matrix_non_20));
<<<<<<< HEAD
norm_session4 = zeros(size(zscore_matrix_non_25));

% Normalize each trial in each session
sessions = {zscore_matrix_non_16, zscore_matrix_non_17, zscore_matrix_non_20, zscore_matrix_non_25};
norm_sessions_nonreward = {norm_session1, norm_session2, norm_session3, norm_session4};
=======

% Normalize each trial in each session
sessions = {zscore_matrix_non_16, zscore_matrix_non_17, zscore_matrix_non_20};
norm_sessions_nonreward = {norm_session1, norm_session2, norm_session3};
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, :);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions_nonreward{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end



<<<<<<< HEAD
%% PLOT - striatum rewarded vs nonrewarded
=======
%% PLOT
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

figure('color','white');
subplot(2, 1, 1)
hold on
ax = gca;
ax.FontSize = 15;

% rewarded
for ii=1:length(norm_sessions_reward)
    current_norm = norm_sessions_reward{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

<<<<<<< HEAD
    cmap = spring(5);
=======
    cmap = spring(4);
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Average Z-score Around Rewards');
end
grid on;
<<<<<<< HEAD
ylim([-1, 2]);
=======
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b
hold off

% nonrewarded
subplot(2, 1, 2)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_nonreward)
    current_norm = norm_sessions_nonreward{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

<<<<<<< HEAD
    cmap = spring(5);
=======
    cmap = spring(4);
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Average Z-score Around non rewarded licks');
end
<<<<<<< HEAD
ylim([-1, 2]);
=======
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b
grid on;
hold off




%% high prob patch vs low prob patch

%% combine sessions

window = 5; % seconds
sampling_rate = 130;
samples = window * sampling_rate;

time = linspace(-window, window, ((samples*2)+1));

% Preallocate matrices for normalized data
norm_session1 = zscore_matrix_high_18;
norm_session2 = zscore_matrix_high_19;
norm_session3 = zscore_matrix_high_22;

norm_sessions_high = {norm_session1, norm_session2, norm_session3};


%% low prob
% Preallocate matrices for normalized data
norm_session1 = zscore_matrix_low_18;
norm_session2 = zscore_matrix_low_19;
norm_session3 = zscore_matrix_low_22;

norm_sessions_low = {norm_session1, norm_session2, norm_session3};

%% PLOT - THIS IS AFTER 30 TRIALS OF EACH PATCH

figure('color','white');
subplot(2, 1, 1)
hold on
ax = gca;
ax.FontSize = 15;
 
% high
for ii=1:length(norm_sessions_high)
    current_norm = norm_sessions_high{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = summer(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Z-score Around Licks in High Probability Patch');
end
grid on;
hold off

% nonrewarded
subplot(2, 1, 2)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_low)
    current_norm = norm_sessions_low{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = summer(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Z-score Around Licks in Low Probability Patch');
end
grid on;
hold off



%% striatum

window = 5; % seconds
sampling_rate = 130;
samples = window * sampling_rate;

time = linspace(-window, window, ((samples*2)+1));

% Preallocate matrices for normalized data
norm_session1 = zscore_matrix_high_16;
norm_session2 = zscore_matrix_high_17;
norm_session3 = zscore_matrix_high_20;

norm_sessions_high = {norm_session1, norm_session2, norm_session3};


%% low prob
% Preallocate matrices for normalized data
norm_session1 = zscore_matrix_low_16;
norm_session2 = zscore_matrix_low_17;
norm_session3 = zscore_matrix_low_20;

norm_sessions_low = {norm_session1, norm_session2, norm_session3};

%% PLOT - THIS IS AFTER 30 TRIALS OF EACH PATCH

figure('color','white');
subplot(2, 1, 1)
hold on
ax = gca;
ax.FontSize = 15;
 
% high
for ii=1:length(norm_sessions_high)
    current_norm = norm_sessions_high{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = spring(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Z-score Around Licks in High Probability Patch');
end
grid on;
hold off

% nonrewarded
subplot(2, 1, 2)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_low)
    current_norm = norm_sessions_low{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = spring(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Z-score Around Licks in Low Probability Patch');
end
grid on;
hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% stay in patch vs switch patch

%% combine sessions

window = 5; % seconds
sampling_rate = 130;
samples = window * sampling_rate;

time = linspace(-window, window, ((samples*2)+1));

% Preallocate matrices for normalized data
norm_session1 = zscore_matrix_stay_18;
norm_session2 = zscore_matrix_stay_19;
norm_session3 = zscore_matrix_stay_22;

norm_sessions_stay = {norm_session1, norm_session2, norm_session3};


%% low prob
% Preallocate matrices for normalized data
norm_session1 = zscore_matrix_switch_18;
norm_session2 = zscore_matrix_switch_19;
norm_session3 = zscore_matrix_switch_22;

norm_sessions_switch = {norm_session1, norm_session2, norm_session3};

%% PLOT 

figure('color','white');
subplot(2, 1, 1)
hold on
ax = gca;
ax.FontSize = 15;
 
% stay
for ii=1:length(norm_sessions_stay)
    current_norm = norm_sessions_stay{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = summer(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Z-score Around Licks When Staying in Patch');
end
grid on;
ylim([-1, 1.5]);
hold off

% switch
subplot(2, 1, 2)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_switch)
    current_norm = norm_sessions_switch{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = summer(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Z-score Around Licks When Switching Patch');
end
grid on;
ylim([-1, 1.5]);
hold off



%% striatum

window = 5; % seconds
sampling_rate = 130;
samples = window * sampling_rate;

time = linspace(-window, window, ((samples*2)+1));

% Preallocate matrices for normalized data
norm_session1 = zscore_matrix_high_16;
norm_session2 = zscore_matrix_high_17;
norm_session3 = zscore_matrix_high_20;

norm_sessions_high = {norm_session1, norm_session2, norm_session3};


%% low prob
% Preallocate matrices for normalized data
norm_session1 = zscore_matrix_low_16;
norm_session2 = zscore_matrix_low_17;
norm_session3 = zscore_matrix_low_20;

norm_sessions_low = {norm_session1, norm_session2, norm_session3};

%% PLOT - THIS IS AFTER 30 TRIALS OF EACH PATCH

figure('color','white');
subplot(2, 1, 1)
hold on
ax = gca;
ax.FontSize = 15;
 
% high
for ii=1:length(norm_sessions_high)
    current_norm = norm_sessions_high{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = spring(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Z-score Around Licks in High Probability Patch');
end
grid on;
hold off

% nonrewarded
subplot(2, 1, 2)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_low)
    current_norm = norm_sessions_low{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 

    cmap = spring(4);
    plot(time, avg_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title('Z-score Around Licks in Low Probability Patch');
end
grid on;
hold off





