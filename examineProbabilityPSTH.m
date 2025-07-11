
%{
function examineProbabilityPSTH

% Look at different types of pokes -
%   1) Rewarded vs unrewarded trials
%   2) Pokes at each of the 7 different ports with high/low reward conditions

sess= {'N7\N7_241216_sess23',...
    }; 

expPath = 'Z:\Buzsakilabspace\LabShare\ZutshiI\patchTask\';

% Prepare PSTH containers
psthReward = cell(13,1);

% Port-to-index mapping (1-3 and 5-7 have high/low reward, port 4 is constant)
portMap = [1, 2, 3, 5, 6, 7, 4];

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath, sess{ii}))
    load(dir('*TrialBehavior.mat').name);
    load(dir('*spikes.cellinfo.mat').name);
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
    load(dir('*cell_metrics.cellinfo.mat').name);

    for kk = 1:length(cell_metrics.UID)
        % Extract pyramidal cells only
        if strcmp(cell_metrics.putativeCellType{kk}, 'Pyramidal Cell')
           
            % Loop through ports and high/low reward conditions
            plotIdx = 1;  
           
            for port = portMap
                for rewardCond = 0:1  % 0 = high reward, 1 = low reward
                   
                    if port == 4 && rewardCond == 1
                        % Skip low reward condition for port 4 (never changes)
                        continue
                    end
                   
                    % Identify trials based on reward condition
                    if rewardCond == 0  % High reward
                        isHighReward = (port <= 3 && behavTrials.patch_number == 0) || ...
                                       (port >= 5 && behavTrials.patch_number == 1);
                    else  % Low reward
                        isHighReward = (port <= 3 && behavTrials.patch_number == 1) || ...
                                       (port >= 5 && behavTrials.patch_number == 0);
                    end
                   
                    % Select trials for this port and reward condition
                    st = behavTrials.timestamps(behavTrials.port == port & isHighReward);
                   
                    % Compute PSTH
                    if ~isempty(st)
                        [stccg, tPSTH] = CCG({spikes.times{kk}, st}, [], 'binSize', 0.1, 'duration', 4, 'norm', 'rate');
                        psthReward{plotIdx} = [psthReward{plotIdx}; stccg(:,2,1)'];
                    else
                        fillArr = nan(1, 41);
                        psthReward{plotIdx} = [psthReward{plotIdx}; fillArr];
                    end
                   
                    plotIdx = plotIdx + 1;
                end
            end
        end
    end
end


%% PLOT

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
YlGnBu = cbrewer('seq', 'YlGnBu', 11);
colormap(YlGnBu)

plotTitles = {
    'Port 1 - High Reward', 'Port 1 - Low Reward', ...
    'Port 2 - High Reward', 'Port 2 - Low Reward', ...
    'Port 3 - High Reward', 'Port 3 - Low Reward', ...
    'Port 4 - Constant Reward', ...
    'Port 5 - High Reward', 'Port 5 - Low Reward', ...
    'Port 6 - High Reward', 'Port 6 - Low Reward', ...
    'Port 7 - High Reward', 'Port 7 - Low Reward'
};

for tt = 1:13
    subplot(4, 4, tt)
   
    if ~isempty(psthReward{tt})
        temp = zscore(psthReward{tt}, [], 2);
        h = imagesc(tPSTH, 1:size(psthReward{tt}, 1), temp);
        set(h, 'AlphaData', ~isnan(temp))
        colorbar
    else
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center')
    end
   
    title(plotTitles{tt})
    xlabel('Time (s)')
    ylabel('Neuron Index')
end

sgtitle('PSTH by Port and Reward Condition')

end

%}


%%%%%%%%%
%{

function examineProbabilityPSTH

% Look at different types of pokes -
%   1) Rewarded vs unrewarded trials
%   2) Pokes at each of the 7 different ports, high vs low reward

sess= {'N7\N7_241216_sess23'};

expPath = 'Z:\Buzsakilabspace\LabShare\ZutshiI\patchTask\';

num_ports = 7;  
psthHigh = cell(1, num_ports);
psthLow = cell(1, num_ports);

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath, sess{ii}))
    file = dir(['*TrialBehavior.mat']);
    load(file(1).name);    
    file = dir(['*spikes.cellinfo.mat']);
    load(file(1).name);      
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);
    
    
    for kk = 1:length(cell_metrics.UID)
        % Determine cell type
        if strcmp(cell_metrics.putativeCellType{kk}, 'Pyramidal Cell') == 1
            cellType = 1;
        else
            cellType = 0;
        end
 
        if cellType == 1  % Only process pyramidal cells
            for trial = 1:behavTrials.num_trials
                patch_number = behavTrials.patch_number(trial);  
                port = behavTrials.port(trial);
                timestamp = behavTrials.timestamps(trial);
    
                % Determine reward status
                if patch_number == 0 && ismember(port, [1, 2, 3])
                    reward_type = 'high';
                elseif patch_number == 1 && ismember(port, [5, 6, 7])
                    reward_type = 'high';
                else
                    reward_type = 'low';
                end
    
                % Skip port 4 (never rewarded)
                if port == 4
                    continue
                end
    
                % Collect spike times for PSTH
                [stccg, tPSTH] = CCG({spikes.times{kk}, timestamp}, [], ...
                    'binSize', 0.1, 'duration', 4, 'norm', 'rate');
                disp(strjoin(string(trial)))
    
                % Store PSTH based on reward type
                if strcmp(reward_type, 'high')
                    psthHigh{port} = [psthHigh{port}; stccg(:, 2, 1)'];
                else
                    psthLow{port} = [psthLow{port}; stccg(:, 2, 1)'];
                end
            end
        end
    end
end

%% PLOT
figure
set(gcf, 'Renderer', 'painters', 'Color', 'w')
YlGnBu = cbrewer('seq', 'YlGnBu', 11);
colormap(YlGnBu)

% Plot high and low reward PSTHs for each port
for port = [1, 2, 3, 5, 6, 7]
    % High reward
    subplot(4, 3, (port-1)*2 - (port > 4)*4 + 1)
    if ~isempty(psthHigh{port})
        temp = zscore(psthHigh{port}, [], 2);
        imagesc(tPSTH, 1:size(psthHigh{port}, 1), temp)
        colorbar
        title(['Port ', num2str(port), ' - High Reward'])
    else
        title(['Port ', num2str(port), ' - No High Reward Trials'])
    end

    % Low reward
    subplot(4, 3, (port-1)*2 - (port > 4)*4 + 2)
    if ~isempty(psthLow{port})
        temp = zscore(psthLow{port}, [], 2);
        imagesc(tPSTH, 1:size(psthLow{port}, 1), temp)
        colorbar
        title(['Port ', num2str(port), ' - Low Reward'])
    else
        title(['Port ', num2str(port), ' - No Low Reward Trials'])
    end
end

% Port 4 (never rewarded)
subplot(4, 3, 13)
title('Port 4 - Never Rewarded')

end


%}


%%%%%%%

function examineProbabilityPSTH

% Look at different types of pokes -
%   1) Rewarded vs unrewarded trials
%   2) Pokes at each of the 7 different ports, high vs low reward

sess= {'N7\N7_241216_sess23'};

expPath = 'Z:\Buzsakilabspace\LabShare\ZutshiI\patchTask\';

num_ports = 7;  
psthHigh = cell(1, num_ports);
psthLow = cell(1, num_ports);

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath, sess{ii}))
    file = dir(['*TrialBehavior.mat']);
    load(file(1).name);    
    file = dir(['*spikes.cellinfo.mat']);
    load(file(1).name);      
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);
    
    
    for kk = 1:length(cell_metrics.UID)
        % Determine cell type
        if strcmp(cell_metrics.putativeCellType{kk}, 'Pyramidal Cell') == 1
            cellType = 1;
        else
            cellType = 0;
        end
 
        if cellType == 1  % Only process pyramidal cells
            for tt = 1:7
                if tt == 1 
                    st_high = behavTrials.timestamps(behavTrials.port == 1 & ...
                        behavTrials.patch_number == 0);
                    st_low = behavTrials.timestamps(behavTrials.port == 1 & ...
                        behavTrials.patch_number == 1);
                elseif tt == 2 
                    st_high = behavTrials.timestamps(behavTrials.port == 2 & ...
                        behavTrials.patch_number == 0);
                    st_low = behavTrials.timestamps(behavTrials.port == 2 & ...
                        behavTrials.patch_number == 1);
                elseif tt == 3 
                    st_high = behavTrials.timestamps(behavTrials.port == 3 & ...
                        behavTrials.patch_number == 0);
                    st_low = behavTrials.timestamps(behavTrials.port == 3 & ...
                        behavTrials.patch_number == 1);
                elseif tt == 4 
                    st = behavTrials.timestamps(behavTrials.port == 4);
                elseif tt == 5 
                    st_high = behavTrials.timestamps(behavTrials.port == 5 & ...
                        behavTrials.patch_number == 1);
                    st_low = behavTrials.timestamps(behavTrials.port == 5 & ...
                        behavTrials.patch_number == 0);
                elseif tt == 6 
                    st_high = behavTrials.timestamps(behavTrials.port == 6 & ...
                        behavTrials.patch_number == 1);
                    st_low = behavTrials.timestamps(behavTrials.port == 6 & ...
                        behavTrials.patch_number == 0);
                elseif tt == 7 
                    st_high = behavTrials.timestamps(behavTrials.port == 7 & ...
                        behavTrials.patch_number == 1);
                    st_low = behavTrials.timestamps(behavTrials.port == 7 & ...
                        behavTrials.patch_number == 0);
                end
    
                if ~isempty(st_high)
                    [stccg, tPSTH_high] = CCG({spikes.times{kk} st_high},[],'binSize',0.1,'duration',4,'norm','rate');
                    psthHigh{tt} = [psthHigh{tt}; stccg(:,2,1)'];               
                else
                    fillArr(1,1:41) = nan;
                    psthHigh{tt} = [psthHigh{tt}; fillArr];               
                end
                if ~isempty(st_low)
                    [stccg, tPSTH_low] = CCG({spikes.times{kk} st_low},[],'binSize',0.1,'duration',4,'norm','rate');
                    psthLow{tt} = [psthLow{tt}; stccg(:,2,1)'];               
                else
                    fillArr(1,1:41) = nan;
                    psthLow{tt} = [psthLow{tt}; fillArr];               
                end
            end
        end
    end
end    

%% PLOT
figure
tiledlayout(2,6);
set(gcf, 'Renderer', 'painters', 'Color', 'w')
YlGnBu = cbrewer('seq', 'YlGnBu', 11);
colormap(YlGnBu)

% Plot high and low reward PSTHs for each port
for port = [1, 2, 3, 5, 6, 7]
    % High reward
    nexttile

    idxT = tPSTH_high<0 & tPSTH_high>=-0.2;
    avgRate1 = nanmean(psthHigh{1}(:,idxT),2);
    avgRate2 = nanmean(psthHigh{2}(:,idxT),2);
    avgRate3 = nanmean(psthHigh{3}(:,idxT),2);
    avgRate5 = nanmean(psthHigh{5}(:,idxT),2);
    avgRate6 = nanmean(psthHigh{6}(:,idxT),2);
    avgRate7 = nanmean(psthHigh{7}(:,idxT),2);
    avgRate = nanmean([avgRate1 avgRate2 avgRate3 avgRate5 avgRate6 avgRate7],2);


    newpsthHigh{port} = psthHigh{port}(avgRate>0.5,:);
    newpsthLow{port} = psthLow{port}(avgRate>0.5,:);
    
    idxT = tPSTH_high<0 & tPSTH_high>=-0.5;
    [~,idxMax2] = max(newpsthHigh{port}(:,idxT),[],2); % 1 for consistent sorting. port for individually sorted
    [~,idxMax] = sort(idxMax2);


    if ~isempty(newpsthHigh{port})
        temp = zscore(newpsthHigh{port}, [], 2);
        h = imagesc(tPSTH_high, 1:size(newpsthHigh{port}, 1), temp(idxMax,:));
        set(h, 'AlphaData', ~isnan(temp))
        colorbar
        title(['Port ', num2str(port), ' - High Reward'])
    else
        title(['Port ', num2str(port), ' - No High Reward Trials'])
    end
end

for port = [1, 2, 3, 5, 6, 7]
    % Low reward
    nexttile
    % comment the following out if you want sorting to stay the same for
    % high/low probability
    
    % idxT = tPSTH_low<0 & tPSTH_low>=-0.2;
    % avgRate1 = nanmean(psthLow{1}(:,idxT),2);
    % avgRate2 = nanmean(psthLow{2}(:,idxT),2);
    % avgRate3 = nanmean(psthLow{3}(:,idxT),2);
    % avgRate5 = nanmean(psthLow{5}(:,idxT),2);
    % avgRate6 = nanmean(psthLow{6}(:,idxT),2);
    % avgRate7 = nanmean(psthLow{7}(:,idxT),2);
    % avgRate = nanmean([avgRate1 avgRate2 avgRate3 avgRate5 avgRate6 avgRate7],2);


    % newpsthLow{port} = psthLow{port}(avgRate>0.5,:);
    % 
    % idxT = tPSTH_low<0 & tPSTH_low>=-0.5;
    % [~,idxMax2] = max(newpsthLow{port}(:,idxT),[],2); % 1 for consistent sorting. port for individually sorted
    % [~,idxMax] = sort(idxMax2);

    if ~isempty(newpsthLow{port})
        temp = zscore(newpsthLow{port}, [], 2);
        h = imagesc(tPSTH_low, 1:size(newpsthLow{port}, 1), temp(idxMax,:));
        set(h, 'AlphaData', ~isnan(temp))
        colorbar
        title(['Port ', num2str(port), ' - Low Reward'])
    else
        title(['Port ', num2str(port), ' - No Low Reward Trials'])
    end
end


% one set of plots
figure('color','white')
tiledlayout(1,6)

for tt = [1, 2, 3, 5, 6, 7]
    % high
    nexttile
    ax = gca;
    ax.FontSize = 15;
    col = [0.047058823529412   0.419607843137255   0.121568627450980];
    meanpsth = mean(newpsthHigh{tt},1);
    stdpsth = std(newpsthHigh{tt},1)./sqrt(size(newpsthHigh{tt},1));         
    hold on
    fill([tPSTH_high; flipud(tPSTH_high)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5); 
    hi = line(tPSTH_high,meanpsth,'LineWidth',1.5,'Color',col);
    
    % low    
    col = [0.654901960784314   0.909803921568627   0.364705882352941];
    meanpsth = mean(newpsthLow{tt},1);
    stdpsth = std(newpsthLow{tt},1)./sqrt(size(newpsthLow{tt},1));         
    fill([tPSTH_low; flipud(tPSTH_low)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5); 
    hi = line(tPSTH_low,meanpsth,'LineWidth',1.5,'Color',col);

    xline(0, '--r', 'LineWidth', 1)
    ylim([1 7])    
    hold off
end 


figure('color','white')
ax = gca;
ax.FontSize = 15;

for tt = 1:6
    pl = tt;
    % high
    subplot(2,6,pl)
    if tt > 3
        tt = tt + 1;
    end
    col = [0.047058823529412   0.419607843137255   0.121568627450980];
    meanpsth = nanmean(newpsthHigh{tt},1);
    stdpsth = nanstd(newpsthHigh{tt},1)./sqrt(size(newpsthHigh{tt},1));         
    hold on
    fill([tPSTH_high; flipud(tPSTH_high)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5); 
    hi = line(tPSTH_high,meanpsth,'LineWidth',1.5,'Color',col);
    xline(0, '--r', 'LineWidth', 1)
    ylim([1 7])
    grid on
    hold off
    
    % low    
    subplot(2,6,pl+6)
    col = [0.654901960784314   0.909803921568627   0.364705882352941];
    meanpsth = nanmean(newpsthLow{tt},1);
    stdpsth = nanstd(newpsthLow{tt},1)./sqrt(size(newpsthLow{tt},1));         
    hold on
    fill([tPSTH_low; flipud(tPSTH_low)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5); 
    hi = line(tPSTH_low,meanpsth,'LineWidth',1.5,'Color',col);
    xline(0, '--r', 'LineWidth', 1)
    grid on
    ylim([1 7])    
    hold off
end 


end


