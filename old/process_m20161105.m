% ferdy 2016-04-27
% 325ms fixation before pre-cue marker
% 25-125ms fixation between pre-cue marker and cue onset
% 100ms cue onset to cue offset
% 500-800ms cue offset to array onset
% 850-1050ms array onset to target dim for long hold trials
% 650-850ms array onset to target dim for long hold trials
% hold shape response window 280-730ms
% release shape response window 280-730ms
% min 280ms before saccade allowed

% evt5 = cue onset
% evt1,2,3,4 = cue offset
% evt6 = array onset and also release shape resp win start
% evt7 = target dim
% evt8 = juice

% 1105 g3-g6, ch65-69 above brain
% 1104 g3-g6, ch65-66 above brain
% 1102 g2-g9, ch65-68 above brain
% 1101 g2-g3, ch65-67 above brain
clear;

sessionRange = 1:4;

sessions = {...
        '20161101', 2, 2, 3;
        '20161102', 2, 8, 4;
        '20161104', 3, 4, 2;
        '20161105', 3, 4, 5
        };
numSessions = size(sessions, 1);

lfpVarName = 'FP098';

numFreqInS400 = 92;
numFreqInS250 = 46;
SPreArrayInRFAll = nan(numSessions, numFreqInS400);
SPreArrayExRFAll = nan(numSessions, numFreqInS400);
SPreDimInRFAll = nan(numSessions, numFreqInS400);
SPreDimExRFAll = nan(numSessions, numFreqInS400);
SPreCueInRFAll = nan(numSessions, numFreqInS250);

nChannels = 32;
startChannel = 65;

supChannelOrig = 6;
deepChannelOrig = 27;

inRFLoc = 3;
exRFLoc = 1;

nLoc = 4;

logRoot = 'gratings_task_eye_tracking_';
cohCueOnsetAll = cell(numSessions, nLoc, nChannels);
cohArrayOnsetAll = cell(numSessions, nLoc, nChannels);
cohArrayOnsetHoldAll = cell(numSessions, nLoc, nChannels);
cohTargetDimAll = cell(numSessions, nLoc, nChannels);

cohLineCueOnsetAll = cell(numSessions, nLoc, nChannels);
cohLineArrayOnsetAll = cell(numSessions, nLoc, nChannels);
cohLineArrayOnsetHoldAll = cell(numSessions, nLoc, nChannels);
cohLineTargetDimAll = cell(numSessions, nLoc, nChannels);

%%
for l = sessionRange
    sessionName = sessions{l,1};
    logStartNum = sessions{l,2};
    numLogsToProcess = sessions{l,3};
    xChannelAlignmentOffset = sessions{l,4};
    fprintf('Processing %s...\n', sessionName);

dataDir = sprintf('C:/Users/Ryan/Documents/MATLAB/gratings-task-data/%s/', ...
        sessionName);
dataFile = sprintf('%s/%s-gmerged-fp.mat', dataDir, sessionName);
load(dataFile);
load('params.mat')

% for the V4 power figure in poster
% deepChannel = deepChannelOrig + xChannelAlignmentOffset;
% lfpVarName = sprintf('FP%03d', startChannel + deepChannel - 1);

fprintf('\tUsing LFP var: %s\n', lfpVarName);

%% process events and sort them into different conditions
usefulEvents = getUsefulEvents(dataDir, logStartNum, numLogsToProcess, ...
        EVT02, EVT04, EVT06, EVT07, EVT14, EVT15, EVT16, EVT08);

% unload these variables to workspace
% (cueOnset, cueOnsetByLoc, ...
%         arrayOnset, arrayOnsetRel, arrayOnsetHold, arrayOnsetByLoc, ...
%         arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
%         arrayOnsetShortHoldByLoc, arrayOnsetLongHoldByLoc, ...
%         targetDim, targetDimByLoc, ...
%         targetDimShortHoldByLoc, targetDimLongHoldByLoc, ...
%         nTrialShortHold, nTrialLongHold, rt);
struct2var(usefulEvents);    

%% LFPs
lfpTs = eval([lfpVarName '_ts']);
lfpTsStep = eval([lfpVarName '_ts_step']);
lfpInd = eval([lfpVarName '_ind']);
Fs = round(1/lfpTsStep);
adjLfp = padNaNsToAdjustLfpOffset(eval(lfpVarName), lfpTs, lfpInd, Fs);
% there are nan's in the data at certain events!!

% lowpass FIR filter - just do the whole thing (slow)
% based on eeglab, use FIR1, with filtorder 3*fix(Fs/hicutoff)
% hiCutoffFreq = 50; % low-pass filter at 50 Hz
% bFirLowPass = fir1(3*fix(params.Fs/hiCutoffFreq), hiCutoffFreq/(params.Fs/2), 'low');

%% process every cell in workspace
allVars = whos;
spikeVarIndices = cellfun(@any, strfind({allVars.name}, 'SPK'));
spikeVars = {allVars(spikeVarIndices).name};

for i = 1:numel(spikeVars)
%     process_spikes_m20161105(sessionName, eval(spikeVars{i}), spikeVars{i}, nLoc, ...
%             arrayOnsetRel, arrayOnsetHold, ...
%             arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
%             targetDim, targetDimByLoc, ...
%             targetDimShortHoldDurByLoc, targetDimLongHoldDurByLoc, ...
%             cueOnset, cueOnsetByLoc);
    % sfc with pulvinar
%     process_sfc_m20161105(sessionName, eval(spikeVars{i}), spikeVars{i}, ...
%             adjLfp, lfpVarName, nLoc, params, ...
%             cueOnset, cueOnsetByLoc, ...
%             arrayOnset, arrayOnsetByLoc, ...
%             arrayOnsetRel, arrayOnsetHold, ...
%             arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
%             arrayOnsetShortHoldByLoc, arrayOnsetLongHoldByLoc, ...
%             targetDim, targetDimByLoc, ...
%             targetDimShortHoldDurByLoc, targetDimLongHoldDurByLoc);
    % local sfc
%     lfpVarName = sprintf('FP%s', spikeVars{i}(4:6));
%     adjLfp = padNaNsToAdjustLfpOffset(eval(lfpVarName), lfpTs, lfpInd, Fs);
%     process_sfc_m20161105(sessionName, eval(spikeVars{i}), spikeVars{i}, ...
%             adjLfp, lfpVarName, nLoc, params, ...
%             cueOnset, cueOnsetByLoc, ...
%             arrayOnset, arrayOnsetByLoc, ...
%             arrayOnsetRel, arrayOnsetHold, ...
%             arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
%             arrayOnsetShortHoldByLoc, arrayOnsetLongHoldByLoc, ...
%             targetDim, targetDimByLoc, ...
%             targetDimShortHoldDurByLoc, targetDimLongHoldDurByLoc);
%     close all;
end

%%
cueOnsetLfpWindow = [1 1.5];
cueOnsetLfpT = -cueOnsetLfpWindow(1) * Fs : cueOnsetLfpWindow(2) * Fs - 1;
lfpAroundCueOnset = createdatamatc(adjLfp, cueOnset, Fs, cueOnsetLfpWindow);
assert(~any(any(isnan(lfpAroundCueOnset))));

arrayOnsetLfpWindow = [1.5 1];
arrayOnsetLfpXLim = [-0.7 0.6] * Fs;
arrayOnsetLfpT = -arrayOnsetLfpWindow(1) * Fs : arrayOnsetLfpWindow(2) * Fs - 1;
lfpAroundArrayOnset = createdatamatc(adjLfp, arrayOnset, Fs, arrayOnsetLfpWindow);
assert(~any(any(isnan(lfpAroundArrayOnset))));

targetDimLfpWindow = [1.5 1];
targetDimLfpXLim = [-0.7 0.4] * Fs;
targetDimLfpT = -targetDimLfpWindow(1) * Fs : targetDimLfpWindow(2) * Fs - 1;
lfpAroundTargetDim = createdatamatc(adjLfp, targetDim, Fs, targetDimLfpWindow);
assert(~any(any(isnan(lfpAroundTargetDim))));

lfpAroundCueOnsetByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetHoldByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetShortHoldByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetShortHoldByLocFilt = cell(nLoc, 1);
lfpAroundArrayOnsetLongHoldByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetLongHoldByLocFilt = cell(nLoc, 1);
lfpAroundArrayOnsetRelByLoc = cell(nLoc, 1);
lfpAroundTargetDimByLoc = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundCueOnsetByLoc{i} = createdatamatc(adjLfp, cueOnsetByLoc{i}, Fs, cueOnsetLfpWindow);
    lfpAroundArrayOnsetHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
    lfpAroundArrayOnsetShortHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetShortHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
%     lfpAroundArrayOnsetShortHoldByLocFilt{i} = filtfilt(bFirLowPass, 1, lfpAroundArrayOnsetShortHoldByLoc{i});
    lfpAroundArrayOnsetLongHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetLongHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
%     lfpAroundArrayOnsetLongHoldByLocFilt{i} = filtfilt(bFirLowPass, 1, lfpAroundArrayOnsetLongHoldByLoc{i});
    lfpAroundArrayOnsetRelByLoc{i} = createdatamatc(adjLfp, arrayOnsetRelByLoc{i}, Fs, arrayOnsetLfpWindow);
    lfpAroundTargetDimByLoc{i} = createdatamatc(adjLfp, targetDimByLoc{i}, Fs, targetDimLfpWindow);
end

% cueOnsetLfpXLim = [-0.4 0.7] * Fs;
% figure;
% hold on;
% for i = 1:size(lfpAroundCueOnset, 2)
%     plot(cueOnsetLfpT, lfpAroundCueOnset(:,i));
% end
% 
% plot(cueOnsetLfpT, mean(lfpAroundCueOnset, 2), 'LineWidth', 5, 'Color', 'k');
% xlim(cueOnsetLfpXLim);
% xlabel('Time from Cue Onset (ms)');
% title(sprintf('Cue Evoked Potential, Raw and Average (N=%d)', size(lfpAroundCueOnset, 2)));

%%
% figure;
% hold on;
% 
% cueOnsetLfpXLim = [-0.4 0.7] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     if isempty(lfpAroundCueOnsetByLoc{i})
%         continue;
%     end
%     fHandles(i) = plot(cueOnsetLfpT, mean(lfpAroundCueOnsetByLoc{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(cueOnsetByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(cueOnsetLfpXLim);
% xlabel('Time from Cue Onset (ms)');
% title(sprintf('Cue Evoked Potential by Location (N=%d)', size(lfpAroundCueOnset, 2)));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');

%%

% arrayOnsetLfpXLim = [-0.7 1] * Fs;
% figure;
% hold on;
% for i = 1:size(lfpAroundArrayOnset, 2)
%     plot(arrayOnsetLfpT, lfpAroundArrayOnset(:,i));
% end
% 
% plot(arrayOnsetLfpT, mean(lfpAroundArrayOnset, 2), 'LineWidth', 5, 'Color', 'k');
% xlim(arrayOnsetLfpXLim);
% xlabel('Time from Array Onset (ms)');
% title(sprintf('Array Evoked Potential (N=%d)', size(lfpAroundArrayOnset, 2)));

%%
% figure;
% hold on;
% 
% arrayOnsetLfpXLim = [-0.7 1] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     if isempty(lfpAroundArrayOnsetHoldByLoc{i})
%         continue;
%     end
%     fHandles(i) = plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetHoldByLoc{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetHoldByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(arrayOnsetLfpXLim);
% xlabel('Time from Array Onset (ms)');
% title(sprintf('Array Evoked Potential by Location -- Hold Target (N=%d)', numel(arrayOnsetHold)));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');

%%
% figure;
% hold on;
% 
% arrayOnsetLfpXLim = [-0.7 1] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     assert(~any(any(isnan(lfpAroundArrayOnsetShortHoldByLoc{i}))));
%     if isempty(lfpAroundArrayOnsetShortHoldByLocFilt{i})
%         continue;
%     end
%     fHandles(i) = plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetShortHoldByLocFilt{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetShortHoldByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(arrayOnsetLfpXLim);
% xlabel('Time from Array Onset (ms)');
% title(sprintf('Array Evoked Potential by Location -- Hold Target Short Dur (N=%d)', nTrialShortHoldDur));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');

%%
% figure;
% hold on;
% 
% arrayOnsetLfpXLim = [-0.7 1] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     assert(~any(any(isnan(lfpAroundArrayOnsetLongHoldByLoc{i}))));
%     if isempty(lfpAroundArrayOnsetLongHoldByLocFilt{i})
%         continue;
%     end
%     fHandles(i) = plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetLongHoldByLocFilt{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetLongHoldByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(arrayOnsetLfpXLim);
% xlabel('Time from Array Onset (ms)');
% title(sprintf('Array Evoked Potential by Location -- Hold Target Long Dur (N=%d)', nTrialLongHoldDur));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');

%%
% figure;
% hold on;
% 
% arrayOnsetLfpXLim = [-0.7 0.6] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     if isempty(lfpAroundArrayOnsetRelByLoc{i})
%         continue;
%     end
%     fHandles(i) = plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetRelByLoc{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetRelByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(arrayOnsetLfpXLim);
% xlabel('Time from Array Onset (ms)');
% title(sprintf('Array Evoked Potential by Location -- Release Target (N=%d)', numel(arrayOnsetRel)));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');
    
%%
% figure;
% hold on;
% for i = 1:size(lfpAroundTargetDim, 2)
%     plot(arrayOnsetLfpT, lfpAroundTargetDim(:,i));
% end
% 
% plot(targetDimLfpT, mean(lfpAroundTargetDim, 2), 'LineWidth', 5, 'Color', 'k');
% xlim(targetDimLfpXLim);
% xlabel('Time from Target Dimming (ms)');
% title(sprintf('Target Dimming Potential (N=%d)', size(lfpAroundTargetDim, 2)));

%%
% figure;
% hold on;
% 
% targetDimLfpXLim = [-0.7 0.4] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     if isempty(lfpAroundTargetDimByLoc{i})
%         continue;
%     end
%     fHandles(i) = plot(targetDimLfpT, mean(lfpAroundTargetDimByLoc{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetHoldByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(targetDimLfpXLim);
% xlabel('Time from Target Dimming (ms)');
% title(sprintf('Target Dimming Evoked Potential by Location (N=%d)', numel(targetDim)));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');

%%
load('params.mat');
params.fpass = [5 50];
params.pad = 2;
specgramSlidingWindow = [0.3 0.05];

%% plot coherence between pul and v4 channels
lfpTs = eval(sprintf('FP%03d_ts', startChannel));
lfpInd = eval(sprintf('FP%03d_ind', startChannel));
lfpTsStep = eval(sprintf('FP%03d_ts_step', startChannel)); % sampling period
v4ChannelData = cell(nChannels, 1);
v4LfpAroundCueOnsetByLoc = cell(nChannels, nLoc);
v4LfpAroundArrayOnsetByLoc = cell(nChannels, nLoc);
v4LfpAroundArrayOnsetHoldByLoc = cell(nChannels, nLoc);
v4LfpAroundTargetDimByLoc = cell(nChannels, nLoc);

cueOnsetLfpWindowNew = [0.15 0.5];
arrayOnsetLfpWindowNew = [0.5 0.15];
targetDimLfpWindowNew = [0.5 0.15];

cueOnsetLfpTNew = -cueOnsetLfpWindowNew(1) * Fs : cueOnsetLfpWindowNew(2) * Fs - 1;
arrayOnsetLfpTNew = -arrayOnsetLfpWindowNew(1) * Fs : arrayOnsetLfpWindowNew(2) * Fs - 1;
targetDimLfpTNew = -targetDimLfpWindowNew(1) * Fs : targetDimLfpWindowNew(2) * Fs - 1;

supChannel = supChannelOrig + xChannelAlignmentOffset;
deepChannel = deepChannelOrig + xChannelAlignmentOffset;

for j = [supChannel deepChannel]%1:nChannels
    channelName1 = sprintf('FP%03d', startChannel + j - 1);
    channelName2 = sprintf('FP%03d', startChannel + j - 1 - 1); % one above
    v4ChannelData{j} = eval(channelName1) - eval(channelName2); % bipolar LFP <----------------
    v4ChannelData{j} = padNaNsToAdjustLfpOffset(v4ChannelData{j}, lfpTs, lfpInd, 1/lfpTsStep);
    for i = [inRFLoc exRFLoc]
        v4LfpAroundCueOnsetByLoc{j,i} = createdatamatc(v4ChannelData{j}, ...
                cueOnsetByLoc{i}, Fs, cueOnsetLfpWindowNew);
        assert(~any(any(isnan(v4LfpAroundCueOnsetByLoc{j,i}))));
        v4LfpAroundArrayOnsetByLoc{j,i} = createdatamatc(v4ChannelData{j}, ...
                arrayOnsetByLoc{i}, Fs, arrayOnsetLfpWindowNew);
        assert(~any(any(isnan(v4LfpAroundArrayOnsetByLoc{j,i}))));
        v4LfpAroundArrayOnsetHoldByLoc{j,i} = createdatamatc(v4ChannelData{j}, ...
                arrayOnsetHoldByLoc{i}, Fs, arrayOnsetLfpWindowNew);
        assert(~any(any(isnan(v4LfpAroundArrayOnsetHoldByLoc{j,i}))));
        v4LfpAroundTargetDimByLoc{j,i} = createdatamatc(v4ChannelData{j}, ...
                targetDimByLoc{i}, Fs, targetDimLfpWindowNew);
        assert(~any(any(isnan(v4LfpAroundTargetDimByLoc{j,i}))));
    end
end
% pulvinar
lfpAroundCueOnsetByLocNew = cell(nLoc, 1);
lfpAroundArrayOnsetByLocNew = cell(nLoc, 1);
lfpAroundArrayOnsetHoldByLocNew = cell(nLoc, 1);
lfpAroundTargetDimByLocNew = cell(nLoc, 1);
for i = [inRFLoc exRFLoc]
    lfpAroundCueOnsetByLocNew{i} = createdatamatc(adjLfp, ...
            cueOnsetByLoc{i}, Fs, cueOnsetLfpWindowNew);
    lfpAroundArrayOnsetByLocNew{i} = createdatamatc(adjLfp, ...
            arrayOnsetByLoc{i}, Fs, arrayOnsetLfpWindowNew);
    lfpAroundArrayOnsetHoldByLocNew{i} = createdatamatc(adjLfp, ...
            arrayOnsetHoldByLoc{i}, Fs, arrayOnsetLfpWindowNew);
    lfpAroundTargetDimByLocNew{i} = createdatamatc(adjLfp, ...
            targetDimByLoc{i}, Fs, targetDimLfpWindowNew);
end

%%
load('params.mat');
params.fpass = [5 50];
params.pad = 2;

cohgramSlidingWindow = [0.3 0.025];
cohLinePreArrayWindowOffset = [-0.2 0];
cohLinePreDimWindowOffset = [-0.2 0];
preArrayLfpSelection = arrayOnsetLfpTNew >= cohLinePreArrayWindowOffset(1)*Fs & ...
        arrayOnsetLfpTNew <= cohLinePreArrayWindowOffset(2)*Fs;
preDimLfpSelection = targetDimLfpTNew >= cohLinePreDimWindowOffset(1)*Fs & ...
        targetDimLfpTNew <= cohLinePreDimWindowOffset(2)*Fs;

for k = [inRFLoc exRFLoc]
    for j = [supChannel deepChannel]%1:nChannels
        fprintf('Processing channel %d...\n', j);
%         [cohCueOnsetAll{l,k,j},phi,S12,~,~,t,f,confC,phistd,Cerr] = cohgramc(...
%                 lfpAroundCueOnsetByLocNew{k}, ...
%                 v4LfpAroundCueOnsetByLoc{j,k}, ...
%                 cohgramSlidingWindow, params);
%         [cohArrayOnsetAll{l,k,j},phi,S12,~,~,t,f,confC,phistd,Cerr] = cohgramc(...
%                 lfpAroundArrayOnsetByLocNew{k}, ...
%                 v4LfpAroundArrayOnsetByLoc{j,k}, ...
%                 cohgramSlidingWindow, params);
%         [cohArrayOnsetHoldAll{l,k,j},phi,S12,~,~,t,f,confC,phistd,Cerr] = cohgramc(...
%                 lfpAroundArrayOnsetHoldByLocNew{k}, ...
%                 v4LfpAroundArrayOnsetHoldByLoc{j,k}, ...
%                 cohgramSlidingWindow, params);
%         [cohTargetDimAll{l,k,j},phi,S12,~,~,t,f,confC,phistd,Cerr] = cohgramc(...
%                 lfpAroundTargetDimByLocNew{k}, ...
%                 v4LfpAroundTargetDimByLoc{j,k}, ...
%                 cohgramSlidingWindow, params);
        % spectrums only
        [cohLineArrayOnsetAll{l,k,j},phi,S12,S1,S2,fSpectrumCohLine] = coherencyc(...
                lfpAroundArrayOnsetByLocNew{k}(preArrayLfpSelection,:), ...
                v4LfpAroundArrayOnsetByLoc{j,k}(preArrayLfpSelection,:), ...
                params);
        [cohLineTargetDimAll{l,k,j},phi,S12,S1,S2,fSpectrumCohLine] = coherencyc(...
                lfpAroundTargetDimByLocNew{k}(preDimLfpSelection,:), ...
                v4LfpAroundTargetDimByLoc{j,k}(preDimLfpSelection,:), ...
                params);
        
%         figure;
%         hold on;
%         plot_matrix(cohAll{k,j}, t - arrayOnsetLfpWindowNew(1), f);
%         plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%         ylim(params.fpass);
%         title(sprintf('Coherence between Pulvinar and V4 Channel %d at P%d', j, k));
    end
end

%% inrf - exrf
% for j = [10 30]
%     figure;
%     hold on;
%     imagesc(t - cueOnsetLfpWindowNew(1), f, (cohCueOnsetAll{l,inRFLoc,j} - cohCueOnsetAll{l,exRFLoc,j})');
%     axis xy;
%     plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%     colorbar;
% %     caxis([-max(abs(caxis)) max(abs(caxis))]);
%     caxis([-0.15 0.15]);
%     colormap(getCoolWarmMap());
%     ylim(params.fpass);
%     title(sprintf('%s Diff Coherence between Pulvinar and V4 Channel %d', sessionName, j));
%     
%     plotFileName = sprintf('%s-PUL-V4-ch%d-diffCoh-postCue.png', sessionName, j);
%     export_fig(plotFileName);
%     
%     figure;
%     hold on;
%     imagesc(t - arrayOnsetLfpWindowNew(1), f, (cohArrayOnsetAll{l,inRFLoc,j} - cohArrayOnsetAll{l,exRFLoc,j})');
%     axis xy;
%     plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%     colorbar;
% %     caxis([-max(abs(caxis)) max(abs(caxis))]);
%     caxis([-0.15 0.15]);
%     colormap(getCoolWarmMap());
%     ylim(params.fpass);
%     title(sprintf('%s Diff Coherence between Pulvinar and V4 Channel %d', sessionName, j));
%     
%     plotFileName = sprintf('%s-PUL-V4-ch%d-diffCoh-preArray.png', sessionName, j);
%     export_fig(plotFileName);
%     
%     figure;
%     hold on;
%     imagesc(t - targetDimLfpWindowNew(1), f, (cohTargetDimAll{l,inRFLoc,j} - cohTargetDimAll{l,exRFLoc,j})');
%     axis xy;
%     plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%     colorbar;
% %     caxis([-max(abs(caxis)) max(abs(caxis))]);
%     caxis([-0.15 0.15]);
%     colormap(getCoolWarmMap());
%     ylim(params.fpass);
%     title(sprintf('%s Diff Coherence between Pulvinar and V4 Channel %d', sessionName, j));
%     
%     plotFileName = sprintf('%s-PUL-V4-ch%d-diffCoh-preDim.png', sessionName, j);
%     export_fig(plotFileName);
%     
%     figure;
%     hold on;
%     plot(fSpectrumCohLine, (cohLineArrayOnsetAll{l,inRFLoc,j} - cohLineArrayOnsetAll{l,exRFLoc,j}));
%     xlim(params.fpass);
%     title(sprintf('%s Diff Coherence between Pulvinar and V4 Channel %d', sessionName, j));
%     
%     plotFileName = sprintf('%s-PUL-V4-ch%d-diffCohLine-preArray.png', sessionName, j);
%     export_fig(plotFileName);
%     
%     figure;
%     hold on;
%     plot(fSpectrumCohLine, (cohLineTargetDimAll{l,inRFLoc,j} - cohLineTargetDimAll{l,exRFLoc,j}));
%     xlim(params.fpass);
%     title(sprintf('%s Diff Coherence between Pulvinar and V4 Channel %d', sessionName, j));
%     
%     plotFileName = sprintf('%s-PUL-V4-ch%d-diffCohLine-preDim.png', sessionName, j);
%     export_fig(plotFileName);
% end


end % end numSessions

%%
% % params.fpass = [30 50];
% % yAxis = [-35 -27];
% figure_tr_inch(8, 6);
% subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
% hold on;
% % plot_vector(mean(SPreCueInRFAll, 1), fCue, 'l', [], 'k');
% h1 = plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1)), ...
%         'Color', [1 0 0], ...
%         'LineWidth', 4);
% % plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1) + std(SPreArrayInRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% % plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1) - std(SPreArrayInRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% h2 = plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1)), ...
%         'Color', [0 0 1], ...
%         'LineWidth', 4);
% % plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1) + std(SPreArrayExRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% % plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1) - std(SPreArrayExRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% h3 = plot(fSpectrum, 10*log10(mean(SPreDimInRFAll, 1)), ...
%         'Color', [1 0 1], ...
%         'LineWidth', 4);
% h4 = plot(fSpectrum, 10*log10(mean(SPreDimExRFAll, 1)), ...
%         'Color', [0 1 1], ...
%         'LineWidth', 4);
% % plot_vector(mean(SPreArrayInRFAll, 1), fSpectrum, 'l', [], 'r');
% % plot_vector(mean(SPreArrayInRFAll, 1) + std(SPreArrayInRFAll, 0, 1) / sqrt(size(SPreArrayInRFAll, 1)), fSpectrum, 'l', [], 'k');
% % plot_vector(mean(SPreArrayInRFAll, 1) - std(SPreArrayInRFAll, 0, 1) / sqrt(size(SPreArrayInRFAll, 1)), fSpectrum, 'l', [], 'k');
% % plot_vector(mean(SPreArrayExRFAll, 1), fSpectrum, 'l', [], 'b');
% % plot_vector(mean(SPreDimInRFAll, 1), fSpectrum, 'l', [], 'm');
% % plot_vector(mean(SPreDimExRFAll, 1), fSpectrum, 'l', [], 'c');
% xlim(params.fpass);
% ylim(yAxis);
% xlabel('Frequency (Hz)');
% ylabel('Power (dB)');
% legend([h1 h2 h3 h4], ...
%         {sprintf('\\color[rgb]{%s}Pre-Stimulus, Attend In RF', num2str(h1.Color))...
%         sprintf('\\color[rgb]{%s}Pre-Stimulus, Attend Out RF', num2str(h2.Color)), ...
%         sprintf('\\color[rgb]{%s}Post-Stimulus, Attend In RF', num2str(h3.Color)), ...
%         sprintf('\\color[rgb]{%s}Post-Stimulus, Attend Out RF', num2str(h4.Color))}, ...
%         'Position', [0.5 0.78 0.4 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
% set(gca, 'FontName', 'Calibri')
% set(gca, 'FontSize', 26);
% set(gca, 'FontWeight', 'bold');
% set(gca, 'Box', 'off');
% set(gca, 'LineWidth', 3);
% set(gcf, 'Color', 'w');
% 
% plotFileName = sprintf('meanSessions-%s-powerComparison-%0.1f-%0.1fHz-n%d.png', ...
%         lfpVarName, params.fpass, numel(sessionRange));
% export_fig(plotFileName);
% 
% stop

%%

% TODO add post-cue 
l = sessionRange(1);
xChannelAlignmentOffset = sessions{l,4};
supChannel = supChannelOrig + xChannelAlignmentOffset;
deepChannel = deepChannelOrig + xChannelAlignmentOffset;

cohPreArrayDiffForAvgSup = nan(numSessions, size(cohArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohPreArrayDiffForAvgDeep = nan(numSessions, size(cohArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohPreDimDiffForAvgSup = nan(numSessions, size(cohTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohTargetDimAll{l,inRFLoc,supChannel}, 2));
cohPreDimDiffForAvgDeep = nan(numSessions, size(cohTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayDiffForAvgSup = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayDiffForAvgDeep = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimDiffForAvgSup = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimDiffForAvgDeep = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));

cohLinePreArrayVsPreDimInRFForAvgSup = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayVsPreDimExRFForAvgSup = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayVsPreDimInRFForAvgDeep = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayVsPreDimExRFForAvgDeep = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));

cohLinePreArrayInRFForAvgSup = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayExRFForAvgSup = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayInRFForAvgDeep = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayExRFForAvgDeep = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimInRFForAvgSup = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimExRFForAvgSup = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimInRFForAvgDeep = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimExRFForAvgDeep = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));

for l = sessionRange
    xChannelAlignmentOffset = sessions{l,4};
    supChannel = supChannelOrig + xChannelAlignmentOffset;
    deepChannel = deepChannelOrig + xChannelAlignmentOffset;

    cohPreArrayDiffForAvgSup(l,:,:) = (cohArrayOnsetAll{l,inRFLoc,supChannel} - cohArrayOnsetAll{l,exRFLoc,supChannel});
    cohPreArrayDiffForAvgDeep(l,:,:) = (cohArrayOnsetAll{l,inRFLoc,deepChannel} - cohArrayOnsetAll{l,exRFLoc,deepChannel});
    cohPreDimDiffForAvgSup(l,:,:) = (cohTargetDimAll{l,inRFLoc,supChannel} - cohTargetDimAll{l,exRFLoc,supChannel});
    cohPreDimDiffForAvgDeep(l,:,:) = (cohTargetDimAll{l,inRFLoc,deepChannel} - cohTargetDimAll{l,exRFLoc,deepChannel});
    
    cohLinePreArrayDiffForAvgSup(l,:,:) = (cohLineArrayOnsetAll{l,inRFLoc,supChannel} - cohLineArrayOnsetAll{l,exRFLoc,supChannel});
    cohLinePreArrayDiffForAvgDeep(l,:,:) = (cohLineArrayOnsetAll{l,inRFLoc,deepChannel} - cohLineArrayOnsetAll{l,exRFLoc,deepChannel});
    cohLinePreDimDiffForAvgSup(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,supChannel} - cohLineTargetDimAll{l,exRFLoc,supChannel});
    cohLinePreDimDiffForAvgDeep(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,deepChannel} - cohLineTargetDimAll{l,exRFLoc,deepChannel});
    
    cohLinePreArrayVsPreDimInRFForAvgSup(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,supChannel} - cohLineArrayOnsetAll{l,inRFLoc,supChannel});
    cohLinePreArrayVsPreDimExRFForAvgSup(l,:,:) = (cohLineTargetDimAll{l,exRFLoc,supChannel} - cohLineArrayOnsetAll{l,exRFLoc,supChannel});
    cohLinePreArrayVsPreDimInRFForAvgDeep(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,deepChannel} - cohLineArrayOnsetAll{l,inRFLoc,deepChannel});
    cohLinePreArrayVsPreDimExRFForAvgDeep(l,:,:) = (cohLineTargetDimAll{l,exRFLoc,deepChannel} - cohLineArrayOnsetAll{l,exRFLoc,deepChannel});
    
    cohLinePreArrayInRFForAvgSup(l,:,:) = (cohLineArrayOnsetAll{l,inRFLoc,supChannel});
    cohLinePreArrayExRFForAvgSup(l,:,:) = (cohLineArrayOnsetAll{l,exRFLoc,supChannel});
    cohLinePreArrayInRFForAvgDeep(l,:,:) = (cohLineArrayOnsetAll{l,inRFLoc,deepChannel});
    cohLinePreArrayExRFForAvgDeep(l,:,:) = (cohLineArrayOnsetAll{l,exRFLoc,deepChannel});
    cohLinePreDimInRFForAvgSup(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,supChannel});
    cohLinePreDimExRFForAvgSup(l,:,:) = (cohLineTargetDimAll{l,exRFLoc,supChannel});
    cohLinePreDimInRFForAvgDeep(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,deepChannel});
    cohLinePreDimExRFForAvgDeep(l,:,:) = (cohLineTargetDimAll{l,exRFLoc,deepChannel});
    
end

meanCohPreArrayDiffSup = squeeze(nanmean(cohPreArrayDiffForAvgSup, 1));
meanCohPreArrayDiffDeep = squeeze(nanmean(cohPreArrayDiffForAvgDeep, 1));
meanCohPreDimDiffSup = squeeze(nanmean(cohPreDimDiffForAvgSup, 1));
meanCohPreDimDiffDeep = squeeze(nanmean(cohPreDimDiffForAvgDeep, 1));

meanCohLinePreArrayDiffSup = squeeze(nanmean(cohLinePreArrayDiffForAvgSup, 1));
meanCohLinePreArrayDiffDeep = squeeze(nanmean(cohLinePreArrayDiffForAvgDeep, 1));
meanCohLinePreDimDiffSup = squeeze(nanmean(cohLinePreDimDiffForAvgSup, 1));
meanCohLinePreDimDiffDeep = squeeze(nanmean(cohLinePreDimDiffForAvgDeep, 1));

meanCohLinePreArrayVsPreDimInRFSup = squeeze(nanmean(cohLinePreArrayVsPreDimInRFForAvgSup, 1));
meanCohLinePreArrayVsPreDimExRFSup = squeeze(nanmean(cohLinePreArrayVsPreDimExRFForAvgSup, 1));
meanCohLinePreArrayVsPreDimInRFDeep = squeeze(nanmean(cohLinePreArrayVsPreDimInRFForAvgDeep, 1));
meanCohLinePreArrayVsPreDimExRFDeep = squeeze(nanmean(cohLinePreArrayVsPreDimExRFForAvgDeep, 1));

meanCohLinePreArrayInRFSup = squeeze(nanmean(cohLinePreArrayInRFForAvgSup));
meanCohLinePreArrayExRFSup = squeeze(nanmean(cohLinePreArrayExRFForAvgSup));
meanCohLinePreArrayInRFDeep = squeeze(nanmean(cohLinePreArrayInRFForAvgDeep));
meanCohLinePreArrayExRFDeep = squeeze(nanmean(cohLinePreArrayExRFForAvgDeep));
meanCohLinePreDimInRFFSup = squeeze(nanmean(cohLinePreDimInRFForAvgSup));
meanCohLinePreDimExRFSup = squeeze(nanmean(cohLinePreDimExRFForAvgSup));
meanCohLinePreDimInRFDeep = squeeze(nanmean(cohLinePreDimInRFForAvgDeep));
meanCohLinePreDimExRFDeep = squeeze(nanmean(cohLinePreDimExRFForAvgDeep));

%%
cohLineYLim = [0 0.4];
cohLineXLim = [5 50];

for l = sessionRange
    sessionName = sessions{l,1};
    xChannelAlignmentOffset = sessions{l,4};
    supChannel = supChannelOrig + xChannelAlignmentOffset;
    deepChannel = deepChannelOrig + xChannelAlignmentOffset;
    
figure_tr_inch(8, 6);
subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, squeeze(cohLinePreArrayInRFForAvgSup(l,:,:)), ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, squeeze(cohLinePreArrayExRFForAvgSup(l,:,:)), ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend({'Attend In RF', ...
        'Attend Out RF'}, ...
        'Position', [0.66 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('%s-%s-V4-ch%d-cohLine-preArray-n%d.png', ...
        sessionName, lfpVarName, supChannel, numel(sessionRange));
export_fig(plotFileName);


figure_tr_inch(8, 6);
subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, squeeze(cohLinePreDimInRFForAvgSup(l,:,:)), ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, squeeze(cohLinePreDimExRFForAvgSup(l,:,:)), ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend({'Attend In RF', ...
        'Attend Out RF'}, ...
        'Position', [0.66 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('%s-%s-V4-ch%d-cohLine-preDim-n%d.png', ...
        sessionName, lfpVarName, supChannel, numel(sessionRange));
export_fig(plotFileName);


figure_tr_inch(8, 6);
subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, squeeze(cohLinePreArrayInRFForAvgDeep(l,:,:)), ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, squeeze(cohLinePreArrayExRFForAvgDeep(l,:,:)), ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend({'Attend In RF', ...
        'Attend Out RF'}, ...
        'Position', [0.66 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('%s-%s-V4-ch%d-cohLine-preArray-n%d.png', ...
        sessionName, lfpVarName, deepChannel, numel(sessionRange));
export_fig(plotFileName);


figure_tr_inch(8, 6);
subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, squeeze(cohLinePreDimInRFForAvgDeep(l,:,:)), ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, squeeze(cohLinePreDimExRFForAvgDeep(l,:,:)), ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend({'Attend In RF', ...
        'Attend Out RF'}, ...
        'Position', [0.66 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('%s-%s-V4-ch%d-cohLine-preDim-n%d.png', ...
        sessionName, lfpVarName, deepChannel, numel(sessionRange));
export_fig(plotFileName);
    
    
end % for each session

%%
cohLineYLim = [0 0.17];

figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, meanCohLinePreArrayInRFSup, ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, meanCohLinePreArrayExRFSup, ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
% xlabel('Frequency (Hz)');
ylabel('Coherence');
% legend({'Attend In RF', ...
%         'Attend Out RF'}, ...
%         'Position', [0.66 0.82 0.3 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preArray-n%d.png', ...
        lfpVarName, supChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
h1 = plot(fSpectrumCohLine, meanCohLinePreDimInRFFSup, ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
h2 = plot(fSpectrumCohLine, meanCohLinePreDimExRFSup, ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
% xlabel('Frequency (Hz)');
% ylabel('Coherence');
legend([h1 h2], ...
        {sprintf('\\color[rgb]{%s}Attend In RF', num2str(h1.Color)), ...
        sprintf('\\color[rgb]{%s}Attend Out RF', num2str(h2.Color))}, ...
        'Position', [0.65 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preDim-n%d.png', ...
        lfpVarName, supChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, meanCohLinePreArrayInRFDeep, ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, meanCohLinePreArrayExRFDeep, ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence');
% legend({'Attend In RF', ...
%         'Attend Out RF'}, ...
%         'Position', [0.66 0.82 0.3 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preArray-n%d.png', ...
        lfpVarName, deepChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, meanCohLinePreDimInRFDeep, ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, meanCohLinePreDimExRFDeep, ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence Difference');
% legend({'Attend In RF', ...
%         'Attend Out RF'}, ...
%         'Position', [0.66 0.82 0.3 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preDim-n%d.png', ...
        lfpVarName, deepChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');

%%
cohLineYLim = [-0.04 0.12];

figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(cohLineXLim, [0 0], 'Color', ones(3,1)*0.1, 'LineWidth', 1);
h1 = plot(fSpectrumCohLine, meanCohLinePreArrayVsPreDimInRFSup, ...
        ':', 'Color', [1 0 0], ...
        'LineWidth', 4);
h2 = plot(fSpectrumCohLine, meanCohLinePreArrayVsPreDimExRFSup, ...
        ':', 'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence Difference');
legend([h1 h2], ...
        {sprintf('\\color[rgb]{%s}Attend In RF (Post-Pre Stimulus)', num2str(h1.Color)), ...
        sprintf('\\color[rgb]{%s}Attend Out RF (Post-Pre Stimulus)', num2str(h2.Color))}, ...
        'Position', [0.45 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preArrayVsPreDim-n%d.png', ...
        lfpVarName, supChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(cohLineXLim, [0 0], 'Color', ones(3,1)*0.1, 'LineWidth', 1);
plot(fSpectrumCohLine, meanCohLinePreArrayVsPreDimInRFDeep, ...
        ':', 'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, meanCohLinePreArrayVsPreDimExRFDeep, ...
        ':', 'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence DIfference');
% legend({'Attend In RF', ...
%         'Attend Out RF'}, ...
%         'Position', [0.66 0.82 0.3 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preArrayVsPreDim-n%d.png', ...
        lfpVarName, deepChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


%%
cohLineYLim = [-0.05 0.07];

figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(cohLineXLim, [0 0], 'Color', ones(3,1)*0.1, 'LineWidth', 1);
plot(fSpectrumCohLine, meanCohLinePreArrayDiffSup, ...
        'Color', [0.6 0.1 0.9], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, meanCohLinePreArrayDiffDeep, ...
        'Color', [0 0.7 0.3], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence Difference');
% legend({'Attend In RF', ...
%         'Attend Out RF'}, ...
%         'Position', [0.66 0.82 0.3 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-chOff%d-diff-cohLine-preArray-n%d.png', ...
        lfpVarName, supChannelOrig, deepChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(cohLineXLim, [0 0], 'Color', ones(3,1)*0.1, 'LineWidth', 1);
h1 = plot(fSpectrumCohLine, meanCohLinePreDimDiffSup, ...
        'Color', [0.6 0.1 0.9], ...
        'LineWidth', 4);
h2 = plot(fSpectrumCohLine, meanCohLinePreDimDiffDeep, ...
        'Color', [0 0.7 0.3], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
% ylabel('Coherence');

set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');
legend([h1 h2], ...
        {sprintf('\\color[rgb]{%s}Pul & Superficial V4', num2str(h1.Color)), ...
        sprintf('\\color[rgb]{%s}Pul & Deep V4', num2str(h2.Color))}, ...
        'Position', [0.6 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-chOff%d-diff-cohLine-preDim-n%d.png', ...
        lfpVarName, supChannelOrig, deepChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


stop

%%
maxAbsCAxis = 0.1;
figure_tr_inch(12,6);

subaxis(1, 2, 1);
hold on;
imagesc(t - arrayOnsetLfpWindowNew(1), f, meanCohPreArrayDiffSup');
axis xy;
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
colorbar;
%     caxis([-max(abs(caxis)) max(abs(caxis))]);
caxis([-maxAbsCAxis maxAbsCAxis]);
colormap(getCoolWarmMap());
ylim(params.fpass);
title(sprintf('Mean Diff Coherence between Pulvinar and V4 Channel %d (n=%d)', ...
        supChannel, numel(sessionRange)));

subaxis(1, 2, 2);
hold on;
imagesc(t - arrayOnsetLfpWindowNew(1), f, meanCohPreArrayDiffDeep');
axis xy;
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
colorbar;
%     caxis([-max(abs(caxis)) max(abs(caxis))]);
caxis([-maxAbsCAxis maxAbsCAxis]);
colormap(getCoolWarmMap());
ylim(params.fpass);
title(sprintf('Mean Diff Coherence between Pulvinar and V4 Channel %d (n=%d)', ...
        deepChannel, numel(sessionRange)));

plotFileName = sprintf('meanSessions-PUL-V4-ch%d-ch%d-meanDiffCoh-preArray-n%d.png', ...
        supChannel, deepChannel, numel(sessionRange));
export_fig(plotFileName);

figure_tr_inch(12,6);

subaxis(1, 2, 1);
hold on;
imagesc(t - arrayOnsetLfpWindowNew(1), f, meanCohPreDimDiffSup');
axis xy;
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
colorbar;
%     caxis([-max(abs(caxis)) max(abs(caxis))]);
caxis([-maxAbsCAxis maxAbsCAxis]);
colormap(getCoolWarmMap());
ylim(params.fpass);
title(sprintf('Mean Diff Coherence between Pulvinar and V4 Channel %d (n=%d)', ...
        supChannel, numel(sessionRange)));

subaxis(1, 2, 2);
hold on;
imagesc(t - arrayOnsetLfpWindowNew(1), f, meanCohPreDimDiffDeep');
axis xy;
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
colorbar;
%     caxis([-max(abs(caxis)) max(abs(caxis))]);
caxis([-maxAbsCAxis maxAbsCAxis]);
colormap(getCoolWarmMap());
ylim(params.fpass);
title(sprintf('Mean Diff Coherence between Pulvinar and V4 Channel %d (n=%d)', ...
        deepChannel, numel(sessionRange)));

plotFileName = sprintf('meanSessions-PUL-V4-ch%d-ch%d-meanDiffCoh-preDim-n%d.png', ...
        supChannel, deepChannel, numel(sessionRange));
export_fig(plotFileName);

