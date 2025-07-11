function [tracking] = getPatchTracking(varargin)

% Gets position tracking for each sub-session and concatenate all of them so they are 
% aligned with LFP and spikes. Default is using the analog input from Bonsai detected position. 
%
%% USAGE
%
%   [tracking] = getToneTracking(varargin)
%
% INPUTS
%   basePath        -(default: pwd) basePath for the recording file, in buzcode format:
%   analogInputPos  - Analog channel that has tracking information.
%   fs              - sampling rate for behavior. default 1250
%   trackImgLength  - Distance of track in the video file, default, 410.
%   trackLength     - Actual length of the track (in cm)
%   freqRange       - Frequency range of the Tone, default 1000 to 22000
%   saveMat         - default true
%   forceReload     - default false
%
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   position.x               - x position in cm/ normalize
%   
%   timestamps      - in seconds, if Basler ttl detected, sync by them

%   HISTORY:
%     - Ipshita Zutshi 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAVE behavTrials from this code, because it will align timestamps

p = inputParser;
addParameter(p,'basepath',pwd,@isstr)
addParameter(p,'fs',150,@isnumeric);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceReload',true,@islogical)

parse(p,varargin{:});
basepath = p.Results.basepath;
fs = p.Results.fs;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;


%% In case tracking already exists 
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) && ~forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

 %% Find subfolder recordings
cd(basepath);
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    count = 1;
    for ii = 1:size(MergePoints.foldernames,2)
         if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'top*']))   
            cd([basepath filesep MergePoints.foldernames{ii}]); 
            fprintf('Computing tracking in %s folder \n',MergePoints.foldernames{ii});
            patchBehav{count} = getPatchBehavior(); 
            behavTracking{count}= PosPatchTracking([],'convFact',0.2732,'forceReload',false); % computing trajectory
            
            % photometry
            photometry_file = dir('*photometry.mat'); %%% might have to improve this method so i make sure i get the behavior photometry
            if ~isempty(photometry_file)
                photom_file_name = load(photometry_file.name); %%%
                photomBehav{count} = getSyncPhotometry_IZ(photom_file_name.photometryData); %%%
                photometry_exists = 1; % 1 if it does, 0 if it doesn't
            else
                photometry_exists = 0;
            end
            
            trackFolder(count) = ii; 
            count = count + 1;
        end
    end
    cd(basepath);
else
    error('missing MergePoints, quiting...');
end

%% Concatenate and sync timestamps
ts = []; subSessions = []; maskSessions = [];
tsBehav= []; % behavior
tsLicks = []; % licks
tsPhot = []; % photometry 
if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    for ii = 1:length(trackFolder)
        if strcmpi(MergePoints.foldernames{trackFolder(ii)},behavTracking{ii}.folder)
            sumTs = behavTracking{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
            ts = [ts; sumTs]; % figure this out

            subSessions = [subSessions; MergePoints.timestamps(trackFolder(ii),1:2)];
            maskSessions = [maskSessions; ones(size(sumTs))*ii];
            
            sumTs = patchBehav{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
            tsBehav = [tsBehav; sumTs];

            if photometry_exists == 1
                sumTsPhot = photomBehav{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1); %%%
                tsPhot = [tsPhot; sumTsPhot]; %%%
            end
            
            %{
            add licks back in if you add licks in behavTrials
            sumTsLick = patchBehav{ii}.lick_timestamps + MergePoints.timestamps(trackFolder(ii),1);
            tsLicks = [tsLicks; sumTsLick];
            %}
            
        else
            error('Folders name does not match!!');
        end
    end
else
    warning('No MergePoints file found. Concatenating timestamps...');
    for ii = 1:length(trackFolder)
        sumTs = max(ts)+ behavTracking{ii}.timestamps;
        subSessions = [subSessions; [sumTs(1) sumTs(end)]];
        ts = [ts; sumTs];
    end
end

% Concatenating tracking fields...
x = []; y = []; vx = []; vy = []; v = []; folder = []; samplingRate = []; description = [];framesDropped = [];
for ii = 1:size(behavTracking,2) 
    x = [x; behavTracking{ii}.position.x];     
    y = [y; behavTracking{ii}.position.y]; 

    vx = [vx; behavTracking{ii}.position.vx];     
    vy = [vy; behavTracking{ii}.position.vy]; 
    v = [v; behavTracking{ii}.position.v];     

    folder{ii} = behavTracking{ii}.folder; 
    samplingRate = [samplingRate; behavTracking{ii}.samplingRate];  
    description{ii} = behavTracking{ii}.description;  
    framesDropped{ii} = behavTracking{ii}.framesDropped;  
end

tracking.position.x = x;
tracking.position.y = y;
tracking.position.vx = vx;
tracking.position.vy = vy;
tracking.position.v = v;

tracking.folders = folder;
tracking.samplingRate = samplingRate;
tracking.timestamps = ts;
tracking.framesDropped = framesDropped;
tracking.events.subSessions =  subSessions;
tracking.events.subSessionsMask = maskSessions;

behavTrials = patchBehav{1};
behavTrials.timestamps = tsBehav;
behavTrials.lick_timestamps = tsLicks;

if photometry_exists == 1
    photometry = photomBehav{1};
    photometry.timestamps = tsPhot;
end

%% save tracking 
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
if saveMat
    save([basepath filesep sessionInfo.FileName '.Tracking.Behavior.mat'],'tracking');
    save([basepath filesep sessionInfo.FileName '.TrialBehavior.mat'],'behavTrials');
    if photometry_exists == 1
        save([basepath filesep sessionInfo.FileName '.PhotometryBehav.mat'],'photometry');
    end
end

end

