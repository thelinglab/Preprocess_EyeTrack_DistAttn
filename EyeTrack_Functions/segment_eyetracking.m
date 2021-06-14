function trial = segment_eyetracking(eyeData,settings)
% segment eye tracker data
%
% Inputs:
% eyeData: struct with relevant eye tracking variables
% settings: settings struct, which specifies segmentation parameters
%
% Outputs:
% trial: structure with segmented data and segmentation parameters
%
% Note: Past versions of this function did not properly handle missing data
% points at the start of the trial. If data points were missing, the first
% available sample was assumed to be the first wanted sample point. So
% instead of there being NaNs at the start of the trial, they ended up at
% the end of the trial. The problem has been corrected in this function.
%
% IMPORTANT FIX: sometimes the recorded time for an event marker
% doesn't line up with the times at which data was sampled. For
% example, when sampling at 500 Hz, we might have data sampled at 1000
% ms, 1002 ms, 1004 ms. However, it's possible that a time-locking
%  message was sent at 1003 ms. In this case, we would try to
% get the time points at 1001 ms, 1003 ms, 1005 ms etc. but these don't
% exist! To deal with this problem, we shift the tWindow index by 1 ms
% if no data was found. This fix isn't necessary for a sampling rate of 1000 Hz. 
% See "fix for when marker time is out of sync with sample points" below

fprintf('segmenting eye tracking data...');
tic

% grab block and trial markers
nMessages = length(eyeData.messages);

blockMessage = {};
blockMessageIdx = [];
trialMessage = {};
trialMessageIdx = [];

% loop through all messages
bCnt = 1;
tCnt = 1;
for m = 1:nMessages

    tmpMessage = eyeData.messages{m};
    if length(tmpMessage) > 5 % all block and trial messages are
        if strcmp('BLOCK',tmpMessage(1:5))
            blockMessage{bCnt} = tmpMessage;
            blockMessageIdx = [blockMessageIdx,m];
            bCnt = bCnt+1;
        end
        if strcmp('TRIAL',tmpMessage(1:5))
            trialMessage{tCnt} = tmpMessage;
            trialMessageIdx = [trialMessageIdx,m];
            tCnt = tCnt+1; 
        end
    end
end

% check that the number of block and trial messages match
if length(blockMessage) ~= length(trialMessage)
    error('the number of block messags and trial messages do not match');
end

trial.nTrials = length(blockMessage); 

% check that time-locking messages are present
timeLockInd = strcmp(settings.seg.timeLockMessage,eyeData.messages); % index the time-locking message (e.g., 'StimOnset')
if sum(timeLockInd) == 0
    error('Dod not find any time-locking events. Did you specify the fight event marker in the settings file?'); 
end

% save times vector for each trial
trial.times = -settings.seg.preTime:eyeData.rateAcq:settings.seg.postTime; % time points in segment
trial.nSamps = length(trial.times); % expected number of samples per segment

% preallocate matrices for segmented data (all the same size)
trial.gx = nan(trial.nTrials,trial.nSamps); trial.gy = trial.gx; trial.pa = trial.gx; trial.exist = trial.gx;
trial.timeLockTimes = nan(trial.nTrials,1); trial.startTimes = trial.timeLockTimes; trial.endTimes = trial.timeLockTimes; 

% create structs to put data in
trial_exist = trial.exist; 
trial_gx = {};
trial_gy = {};
trial_pa = {};

% loop through trials and segment data
for t = 1:trial.nTrials
    
     trial_gx = []; trial_gy = []; trial_pa = []; existInd = [];
    
    % grab vector of messages for the curent trial
    if t < trial.nTrials
        trialIdx = blockMessageIdx(t):blockMessageIdx(t+1);
    else
        trialIdx = blockMessageIdx(t):nMessages;
    end
    
    markTrial = zeros(1,nMessages);
    markTrial(trialIdx) = 1;
    
    lockMessageIdx = zeros(1,nMessages);
    lockMessageIdx(markTrial == 1 & timeLockInd == 1) = 1;
    lockMessageIdx = logical(lockMessageIdx);
        
    % if there is a lock message for the current trial, grab data
    if sum(lockMessageIdx) > 0
        
        % specify start and end times of segment
        trial.timeLockTimes(t) = eyeData.eventTimes(lockMessageIdx); % times for time-locking messsage
        trial.startTimes(t) = double(trial.timeLockTimes(t)) - settings.seg.preTime; % start time, ms
        trial.endTimes(t) = double(trial.timeLockTimes(t)) + settings.seg.postTime;  % end time, ms
        
        % grab the start and end of trial t
        tStart = trial.startTimes(t); tEnd = trial.endTimes(t);
        
        % specify window of interest
        tWindow = tStart:double(eyeData.rateAcq):tEnd;
        
        % index times of interest with logical
        tWindowInd = ismember(double(eyeData.sampleTimes),tWindow);
        
        % fix for when marker time is out of sync with sample points
        if eyeData.rateAcq == 2
            if sum(tWindowInd) == 0
                tWindow = tWindow-1;
                tWindowInd = ismember(double(eyeData.sampleTimes),tWindow);
            end
        end
        
        % throw an error if sampling rate is less than 500 Hz
        if eyeData.rateAcq > 2
            error('Sampling rate lower than 500 Hz. Have not prepared the fix above for sampling freqs lower than 500 Hz')
        end
        
        % create index of the time points that actually exist in the data (i.e., that were recorded).
        existInd = ismember(tWindow,double(eyeData.sampleTimes)); % FIXED - changed to double
        
        % determine which eye was recorded for trial t
        recordedEye = eyeData.RecordedEye(t);
        
        % grab the relevant segment of data (from the recorded eye)
        trial_gx = eyeData.gx(recordedEye,tWindowInd);
        trial_gy = eyeData.gy(recordedEye,tWindowInd);
        trial_pa = eyeData.pa(recordedEye,tWindowInd);
        
        % save exist to the trial structure to make it easy to check where data is missing
        trial.exist(t,:) = existInd;
        
        % save to matrix;
        trial.gx(t,existInd) = trial_gx;
        trial.gy(t,existInd) = trial_gy;
        trial.pa(t,existInd) = trial_pa;
        
    end
    
end

% save vector of block and trials numbers
for t = 1:trial.nTrials
    
    tmp = blockMessage{t};
    tmp = tmp(7:end);
    trial.blockNum(t) = str2num(tmp); 
    
    tmp = trialMessage{t};
    tmp = tmp(7:end);
    trial.trialNum(t) = str2num(tmp); 
    
end

% plot the missing data to alert experimenter to problems
figure; imagesc(trial.exist);
title('Missing samples (1 = data present, 0 = data missing)')
xlabel('Samples')
ylabel('Trials')
caxis([0 1]) % FIXED - hard code color limit
colorbar

fprintf('\n')
toc