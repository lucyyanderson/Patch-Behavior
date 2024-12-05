function [patch] = getPatchBehavior(varargin)

p = inputParser;
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'plotfig',true,@islogical)
addParameter(p,'forceRun',true,@islogical)
addParameter(p,'updatedIntan',true,@islogical)

parse(p,varargin{:});
saveMat = p.Results.saveMat;
plotfig = p.Results.plotfig;
forceRun = p.Results.forceRun;
updatedIntan = p.Results.updatedIntan;

basepath = pwd;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.TrialBehavior.Events.mat'])) && ~forceRun
    disp('Trial behavior already detected! Loading file.');
    file = dir([basepath filesep '*.TrialBehavior.Events.mat']);
    load(file.name);
    return
end

%% Get digital inputs
if exist('settings.xml')
    delete 'settings.xml'
end

disp('Loading digital In...');
digitalIn = bz_getDigitalIn; 
if isempty(digitalIn)
    patchBehav = [];
    return
end

% add columns to each data with their port number
lick_1 = [digitalIn.timestampsOn{1,3}, ones(length(digitalIn.timestampsOn{1,3}),1)];
lick_2 = [digitalIn.timestampsOn{1,4}, 2*ones(length(digitalIn.timestampsOn{1,4}),1)];
lick_3 = [digitalIn.timestampsOn{1,5}, 3*ones(length(digitalIn.timestampsOn{1,5}),1)];
lick_4 = [digitalIn.timestampsOn{1,6}, 4*ones(length(digitalIn.timestampsOn{1,6}),1)];
lick_5 = [digitalIn.timestampsOn{1,7}, 5*ones(length(digitalIn.timestampsOn{1,7}),1)];
lick_6 = [digitalIn.timestampsOn{1,8}, 6*ones(length(digitalIn.timestampsOn{1,8}),1)];
lick_7 = [digitalIn.timestampsOn{1,9}, 7*ones(length(digitalIn.timestampsOn{1,9}),1)];

% creat matrix from individual arrays
all_starts = [lick_1; lick_2; lick_3; lick_4; lick_5; lick_6; lick_7];
all_starts = sortrows(all_starts);    

trial_licks = []; % array which stores the time of the lick indicating a trial. so it excludes multiple licks at one port, and only counts the first lick (which is the one that would be rewarded)
choice_port = [];
most_recent_lick = 100;

% i could add a num_licks to this for loop to create an array that stores
% how many times the mouse licks at each port. not sure if it's totally
% necessary to do rn, but i could if we want that info
for j = 1:length(all_starts)
    start_time = all_starts(j, 1); % gets nonzero number in the row
    current_port = all_starts(j, 2); 
    if most_recent_lick == current_port
        continue
    else
        most_recent_lick = current_port;
        trial_licks = [trial_licks; start_time]; 
        choice_port = [choice_port; current_port];
    end
end

%% Initialize
patch.timestamps = trial_licks;
patch.port =  choice_port;
patch.reward_outcome = zeros(size(trial_licks,1),1); % rewarded (1), not rewarded (0)
patch.patch_probability = zeros(size(trial_licks,1),7); % probabilities of all ports
patch.patch_numbers = zeros(size(trial_licks,1),1); % what patch has higher rewards
patch.patch_trials = zeros(size(trial_licks,1),1); % what number trial it is within the patch

%% Get solenoid data
sol = 0;
for i = 1:length(trial_licks)
    time = trial_licks(i);
    prt = choice_port(i);
    if prt == 1
        sol = 10;
    elseif prt == 2
        sol = 11;
    elseif prt == 3
        sol = 12;
    elseif prt == 5 % skipping middle port
        sol = 14;
    elseif prt == 6
        sol = 15;
    elseif prt == 7
        sol = 16;
    end
    vals = (digitalIn.timestampsOn{1,sol} > (time - 0.01)) & (digitalIn.timestampsOn{1,sol} < (time + 0.01));
    if vals == 1
        patch.reward_outcome(i) = 1;
    else
        patch.reward_outcome(i) = 0;
    end
end


%% Analyze text file
% this only gets info on the probabilities and patch, because reward and lick data is
% coming from intan. to read all info from the text file, use
% analyze_patch_data.m

filename = 'N5_10.22.txt'; % change each time you analyze a new session
fid = fopen(filename, 'r');

% Initialize variables
trial_count = 0;
trial_numbers = [];

while ~feof(fid)
    % Read first line (trial data)
    trial_line = fgetl(fid);
    trial_count = trial_count + 1;  

    trial_info = split(trial_line, '_');
    trial_numbers(trial_count) = str2double(trial_info{2});  
    current_rewarded_count = str2double(trial_info{3});     
    
    % get reward probabilities
    prob_str = trial_info{5};
    reward_probs = (split(prob_str, ','));
    patch.patch_probability(trial_count, 1) = str2double(reward_probs{1});
    patch.patch_probability(trial_count, 2) = str2double(reward_probs{2});
    patch.patch_probability(trial_count, 3) = str2double(reward_probs{3});
    patch.patch_probability(trial_count, 4) = 0;
    patch.patch_probability(trial_count, 5) = str2double(reward_probs{4});
    patch.patch_probability(trial_count, 6) = str2double(reward_probs{5});
    patch.patch_probability(trial_count, 7) = str2double(reward_probs{6});

    
    % read next line (patch info)

    patch_line = fgetl(fid);
    patch_info = split(patch_line, '_');
    patch.patch_trials(trial_count) = str2double(patch_info{1});
    patch.patch_numbers(trial_count) = str2double(patch_info{2});
        
end

fclose(fid);


%% Plot

lick_matrix = zeros(num_trials,7);

% fill matrix
for k = 1:num_trials
    pt = behavTrials.port(k);
    if behavTrials.reward_outcome(k) == 1
        lick_matrix(k, pt) = 2; % rewarded
    elseif behavTrials.reward_outcome(k) == 0
        lick_matrix(k, pt) = 1; % not rewarded
    end
    % unlicked ports remain 0
end
behavTrials.lick_matrix = lick_matrix;


patch_licks = zeros(3,1);
patch_licks(1) = nnz(behavTrials.port == 1) + nnz(behavTrials.port == 2) + nnz(behavTrials.port == 3);
patch_licks(2) = nnz(behavTrials.port == 4);
patch_licks(3) = nnz(behavTrials.port == 5) + nnz(behavTrials.port == 6) + nnz(behavTrials.port == 7);


figure
subplot(1,3,[1,2])
ylabel('Trial #');
xlabel('Reward Port');
heatmap(lick_matrix);

subplot(1,3,3)
title('Distribution of Licks across patches');
x = ["Patch 0" "Middle port" "Patch 1"];
bar(x, patch_licks);

if saveMat
    C = strsplit(pwd,'\');
    save([basepath filesep C{end} '.TrialBehavior.Behavior.mat'],'behavTrials');
end

disp('done!');

%% Plot how often they lick in the correct patch

lick_outcome = zeros(3, 1);
% first row of lick outcome = rewarded patch, second is middle port,
% third is unrewarded patch

for i = 1:length(patch.patch_number)
    % determine which patch has higher rewards
    if patch.patch_number(i) == 0
        rewarded_ports = [1, 2, 3];
        nonrewarded_ports = [5, 6, 7];
    else 
        rewarded_ports = [5, 6, 7];
        nonrewarded_ports = [1, 2, 3];
    end

    % determine if mouse was licking high or low patch
    if ismember(patch.port(i), rewarded_ports)
        lick_outcome(1) = lick_outcome(1)+1;
    elseif ismember(patch.port(i), nonrewarded_ports)
        lick_outcome(3) = lick_outcome(3)+1;
    else 
        lick_outcome(2) = lick_outcome(2)+1;
    end
end

figure(2)
title('Distribution of Licks across patches');
z = ["Rewarded Patch" "Middle port" "Nonrewarded Patch"];
bar(z, lick_outcome);





