% evt5 = cue onset
% evt1,2,3,4 = cue offset
% evt6 = array onset and also release shape resp win start
% evt7 = target dim
% evt8 = juice

event5 = EVT05(:,1);
event6 = EVT06(:,1);
event7 = EVT07(:,1);
event8 = EVT08(:,1);

% get rid of extra juice pulse at end of block
event8(188) = [];

% ferdy 2016-04-27
% 500-600ms fixation before cue onset
% 100ms cue onset to cue offset
% 500-800ms cue offset to array onset
% 600-1400ms array onset to target dim
% hold shape response window 200-800ms
% release shape response window 200-800ms
% min 280ms before saccade allowed

%%
trialParams = readPresentationLogCorrectTrial('C:\Users\Ryan\Documents\MATLAB\gratings-task-data\f20160427\', 'gratings_task_eye_tracking_', 1);
trialParamsRT = trialParams(:,7);

trialsRecorded = 297; % missed the first some
trialParams = trialParams(end-trialsRecorded:end,:);
assert(size(trialParams,1) == numel(event8));

% get only trials that are not repeats
notRepeatLogical = trialParams(:,3) == 1;
event8 = event8(notRepeatLogical,:);
trialParams = trialParams(notRepeatLogical,:);

trialCueLoc = (trialParams(:,4) + 3) / 2;
% here, P1 is btm right, P2 is top right, P3 is top left, P4 is btm left
holdDur = trialParams(:,6);

%%
% get EVT6 corresponding to correct trials only
maxArrayOnsetToJuiceTime = 2.5; % seconds
maxEvent6ToEvent6Time = 0.3;
arrayOnset = nan(numel(event8), 1);
for i = 1:numel(event8)  
    prevEvent6Ind = find((event8(i) - event6 < maxArrayOnsetToJuiceTime) & (event8(i) - event6 > 0), 1, 'last');
    if ~isempty(prevEvent6Ind)
        arrayOnset(i) = event6(prevEvent6Ind);
    else
        warning('unpaired event8: %d, %f', i, event8(i));
    end
    
    % there are two event6's for each release shape array onset, separated 
    % by about 200ms. use the first one. the other marks the start of the
    % response window.
    event6BeforeEvent6Ind = find((event6(prevEvent6Ind) - event6 < maxEvent6ToEvent6Time) & (event6(prevEvent6Ind) - event6 > 0), 1, 'first');
    if ~isempty(event6BeforeEvent6Ind)
        arrayOnset(i) = event6(event6BeforeEvent6Ind);
    end
end

% split release and hold shapes - not the best way, but works for now
maxArrayOnsetToJuiceTimeReleaseShape = 0.7;
isHoldTrial = event8 - arrayOnset >= maxArrayOnsetToJuiceTimeReleaseShape;
arrayOnsetRel = arrayOnset(~isHoldTrial);
arrayOnsetHold = arrayOnset(isHoldTrial);

nLoc = 4;
arrayOnsetRelByLoc = cell(nLoc,1);
arrayOnsetHoldByLoc = cell(nLoc,1);
holdDurMid = (600+1400)/2;
nTrialShortHoldDur = sum(isHoldTrial & holdDur < holdDurMid);
nTrialLongHoldDur = sum(isHoldTrial & holdDur >= holdDurMid);
arrayOnsetShortHoldByLoc = cell(nLoc,1);
arrayOnsetLongHoldByLoc = cell(nLoc,1);
for i = 1:nLoc
    arrayOnsetRelByLoc{i} = arrayOnset(~isHoldTrial & trialCueLoc == i);
    arrayOnsetHoldByLoc{i} = arrayOnset(isHoldTrial & trialCueLoc == i);
    arrayOnsetShortHoldByLoc{i} = arrayOnset(isHoldTrial & trialCueLoc == i & holdDur < holdDurMid);
    arrayOnsetLongHoldByLoc{i} = arrayOnset(isHoldTrial & trialCueLoc == i & holdDur >= holdDurMid);
end

assert(all(isHoldTrial == (trialParams(:,6) ~= -1)));



%% align spikes to events
spikeVar = SPKC097a;
periArrayOnsetWindow = [1 1]; % seconds before, after

arrayOnsetRelSpikeTimes = createnonemptydatamatpt(spikeVar, arrayOnsetRel, periArrayOnsetWindow);
arrayOnsetHoldSpikeTimes = createnonemptydatamatpt(spikeVar, arrayOnsetHold, periArrayOnsetWindow);

arrayOnsetRelSpikeTimesByLoc = cell(nLoc,1);
arrayOnsetHoldSpikeTimesByLoc = cell(nLoc,1);
for i = 1:nLoc
    arrayOnsetRelSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, arrayOnsetRelByLoc{i}, periArrayOnsetWindow);
    arrayOnsetHoldSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, arrayOnsetHoldByLoc{i}, periArrayOnsetWindow);
end

%% set ylim
origYLim = [0 22];

%% spdf
kernelSigma = 0.02;
arrayOnsetSpdfWindowOffset = [-0.3 1];
arrayOnsetT = computeTForSpdf(periArrayOnsetWindow(1), arrayOnsetSpdfWindowOffset, kernelSigma);

arrayOnsetRelSpdf = edpsth_notranspose(arrayOnsetRelSpikeTimes, kernelSigma, 'n', [], 0, arrayOnsetT);
arrayOnsetHoldSpdf = edpsth_notranspose(arrayOnsetHoldSpikeTimes, kernelSigma, 'n', [], 0, arrayOnsetT);

figure;
hold on;
plot(arrayOnsetT - periArrayOnsetWindow(1), arrayOnsetHoldSpdf);
plot(arrayOnsetT - periArrayOnsetWindow(1), arrayOnsetRelSpdf);
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim([-0.3 0.4]);
xlabel('Time from Array Onset (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('f20160427 - Response to Array Onset (N=%d)', numel(event8)));
legend({sprintf('Target: Hold (N=%d)', numel(arrayOnsetHold)), ...
        sprintf('Target: Release (N=%d)', numel(arrayOnsetRel))}, ...
        'Location', 'NorthWest');

%%
figure;
hold on;

arrayOnsetSpdfWindowOffset = [-0.3 1];
arrayOnsetHoldSpdfByLoc = cell(nLoc,1);
legendEntry = cell(4,1);
for i = 1:nLoc
    arrayOnsetHoldSpdfByLoc{i} = edpsth_notranspose(arrayOnsetHoldSpikeTimesByLoc{i}, kernelSigma, 'n', [], 0, arrayOnsetT);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetHoldByLoc{i}));
    plot(arrayOnsetT - periArrayOnsetWindow(1), arrayOnsetHoldSpdfByLoc{i});
end

plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(arrayOnsetSpdfWindowOffset);
xlabel('Time from Array Onset (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('f20160427 - Response to Array Onset - Hold Shape by Location (N=%d)', numel(arrayOnsetHold)));
legend(legendEntry, ...
        'Location', 'NorthWest');
    
%%
figure;
hold on;

arrayOnsetSpdfWindowOffset = [-0.3 0.7];
arrayOnsetRelSpdfByLoc = cell(nLoc,1);
legendEntry = cell(4,1);
for i = 1:nLoc
    arrayOnsetRelSpdfByLoc{i} = edpsth_notranspose(arrayOnsetRelSpikeTimesByLoc{i}, kernelSigma, 'n', [], 0, arrayOnsetT);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetRelByLoc{i}));
    plot(arrayOnsetT - periArrayOnsetWindow(1), arrayOnsetRelSpdfByLoc{i});
end

plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(arrayOnsetSpdfWindowOffset);
xlabel('Time from Array Onset (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('f20160427 - Response to Array Onset - Release Shape by Location (N=%d)', numel(arrayOnsetRel)));
legend(legendEntry, ...
        'Location', 'NorthWest');

%%
% get EVT7 corresponding to correct trials only
maxTargetDimToJuiceTime = 1; % seconds
targetDimMatch = nan(numel(event8), 1);
rt = nan(numel(event8), 1);
for i = 1:numel(event8)  
    prevEvent7Ind = find((event8(i) - event7 < maxTargetDimToJuiceTime) & (event8(i) - event7 > 0), 1, 'last');
    if ~isempty(prevEvent7Ind)
        targetDimMatch(i) = event7(prevEvent7Ind);
        rt(i) = event8(i) - targetDimMatch(i);
    else
        rt(i) = event8(i) - arrayOnset(i);
    end
end

targetDim = targetDimMatch(~isnan(targetDimMatch));
assert(numel(unique(targetDim)) == numel(targetDim));

targetDimByLoc = cell(nLoc,1);
targetDimShortHoldDurByLoc = cell(nLoc,1);
targetDimLongHoldDurByLoc = cell(nLoc,1);
nTrialShortHoldDur = sum(~isnan(targetDimMatch) & holdDur < holdDurMid);
nTrialLongHoldDur = sum(~isnan(targetDimMatch) & holdDur >= holdDurMid);
for i = 1:nLoc
    targetDimByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & trialCueLoc == i);
    targetDimShortHoldDurByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & trialCueLoc == i & holdDur < holdDurMid);
    targetDimLongHoldDurByLoc{i} = targetDimMatch(~isnan(targetDimMatch) & trialCueLoc == i  & holdDur >= holdDurMid);
end

%% align spikes to events
periTargetDimWindow = [1 1];
targetDimSpikeTimes = createnonemptydatamatpt(spikeVar, targetDim, periTargetDimWindow);

targetDimSpikeTimesByLoc = cell(nLoc,1);
targetDimShortHoldDurSpikeTimesByLoc = cell(nLoc,1);
targetDimLongHoldDurSpikeTimesByLoc = cell(nLoc,1);
for i = 1:nLoc
    targetDimSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimByLoc{i}, periTargetDimWindow);
    targetDimShortHoldDurSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimShortHoldDurByLoc{i}, periTargetDimWindow);
    targetDimLongHoldDurSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimLongHoldDurByLoc{i}, periTargetDimWindow);
end

%% spdf
targetDimSpdfWindowOffset = [-1 0.7];
targetDimT = computeTForSpdf(periTargetDimWindow(1), targetDimSpdfWindowOffset, kernelSigma);

targetDimSpdf = edpsth_notranspose(targetDimSpikeTimes, kernelSigma, 'n', [], 0, targetDimT);

figure;
hold on;
plot(targetDimT - periTargetDimWindow(1), targetDimSpdf);
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('f20160427 - Response to Hold Shape Target Dimming (N=%d)', numel(targetDim)));

%%
figure;
hold on;

targetDimSpdfByLoc = cell(nLoc,1);
legendEntry = cell(4,1);
for i = 1:nLoc
    targetDimSpdfByLoc{i} = edpsth_notranspose(targetDimSpikeTimesByLoc{i}, kernelSigma, 'n', [], 0, targetDimT);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(targetDimByLoc{i}));
    plot(targetDimT - periTargetDimWindow(1), targetDimSpdfByLoc{i});
end

plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('f20160427 - Response to Hold Shape Dimming by Location (N=%d)', numel(targetDim)));
legend(legendEntry, ...
        'Location', 'NorthWest');

%%
figure;
hold on;

targetDimShortHoldDurSpdfByLoc = cell(nLoc,1);
legendEntry = cell(4,1);
for i = 1:nLoc
    targetDimShortHoldDurSpdfByLoc{i} = edpsth_notranspose(targetDimShortHoldDurSpikeTimesByLoc{i}, kernelSigma, 'n', [], 0, targetDimT);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(targetDimShortHoldDurByLoc{i}));
    plot(targetDimT - periTargetDimWindow(1), targetDimShortHoldDurSpdfByLoc{i});
end

plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('f20160427 - Response to Hold Shape Dimming by Location SHORT Hold (N=%d)', nTrialShortHoldDur));
legend(legendEntry, ...
        'Location', 'NorthWest');

%%
figure;
hold on;

targetDimLongHoldDurSpdfByLoc = cell(nLoc,1);
legendEntry = cell(4,1);
for i = 1:nLoc
    targetDimLongHoldDurSpdfByLoc{i} = edpsth_notranspose(targetDimLongHoldDurSpikeTimesByLoc{i}, kernelSigma, 'n', [], 0, targetDimT);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(targetDimLongHoldDurByLoc{i}));
    plot(targetDimT - periTargetDimWindow(1), targetDimLongHoldDurSpdfByLoc{i});
end

plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('f20160427 - Response to Hold Shape Dimming by Location LONG Hold (N=%d)', nTrialLongHoldDur));
legend(legendEntry, ...
        'Location', 'NorthWest');

%%
% get EVT5 corresponding to correct trials only
maxCueOnsetToJuiceTime = 5; % seconds
cueOnset = nan(numel(event8), 1);
for i = 1:numel(event8)  
    prevEvent5Ind = find((event8(i) - event5 < maxCueOnsetToJuiceTime) & (event8(i) - event5 > 0), 1, 'last');
    if ~isempty(prevEvent5Ind)
        cueOnset(i) = event5(prevEvent5Ind);
    else
        warning('unpaired event8: %d, %f', i, event8(i));
    end
end

cueOnsetByLoc = cell(nLoc,1);
for i = 1:nLoc
    cueOnsetByLoc{i} = cueOnset(trialCueLoc == i);
end

%% align spikes to events
periCueOnsetWindow = [1 1];
cueOnsetSpikeTimes = createnonemptydatamatpt(spikeVar, cueOnset, periCueOnsetWindow);

cueOnsetSpikeTimesByLoc = cell(nLoc,1);
for i = 1:nLoc
    cueOnsetSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, cueOnsetByLoc{i}, periCueOnsetWindow);
end

%% spdf
cueOnsetSpdfWindowOffset = [-0.4 0.5];
cueOnsetT = computeTForSpdf(periCueOnsetWindow(1), cueOnsetSpdfWindowOffset, kernelSigma);

cueOnsetSpdf = edpsth_notranspose(cueOnsetSpikeTimes, kernelSigma, 'n', [], 0, cueOnsetT);

figure;
hold on;
plot(cueOnsetT - periCueOnsetWindow(1), cueOnsetSpdf);
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(cueOnsetSpdfWindowOffset);
xlabel('Time from Cue Onset (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('f20160427 - Response to Cue Onset (all locations) (N=%d)', numel(cueOnset)));

%%
figure;
hold on;

cueOnsetSpdfByLoc = cell(nLoc,1);
legendEntry = cell(4,1);
for i = 1:nLoc
    cueOnsetSpdfByLoc{i} = edpsth_notranspose(cueOnsetSpikeTimesByLoc{i}, kernelSigma, 'n', [], 0, cueOnsetT);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(cueOnsetByLoc{i}));
    plot(cueOnsetT - periCueOnsetWindow(1), cueOnsetSpdfByLoc{i});
end

plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(cueOnsetSpdfWindowOffset);
xlabel('Time from Cue Onset (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('f20160427 - Response to Cue Onset by Location (N=%d)', numel(cueOnset)));
legend(legendEntry, ...
        'Location', 'NorthWest');

%%
Fs = round(1/FP097_ts_step);
adjLfp = padNaNsToAdjustLfpOffset(FP097, FP097_ts, FP097_ind, Fs);
% there are nan's in the data at certain events!!

% lowpass FIR filter - just do the whole thing (slow)
% based on eeglab, use FIR1, with filtorder 3*fix(Fs/hicutoff)
hiCutoffFreq = 50; % low-pass filter at 50 Hz
bFirLowPass = fir1(3*fix(params.Fs/hiCutoffFreq), hiCutoffFreq/(params.Fs/2), 'low');

%%
cueOnsetLfpWindow = [1 1.5];
cueOnsetLfpXLim = [-0.4 0.7] * Fs;
cueOnsetLfpT = -cueOnsetLfpWindow(1) * Fs : cueOnsetLfpWindow(2) * Fs - 1;
lfpAroundCueOnset = createdatamatc(adjLfp, cueOnset, Fs, cueOnsetLfpWindow);
assert(~any(any(isnan(lfpAroundCueOnset))));

figure;
hold on;
for i = 1:size(lfpAroundCueOnset, 2)
    plot(cueOnsetLfpT, lfpAroundCueOnset(:,i));
end

plot(cueOnsetLfpT, mean(lfpAroundCueOnset, 2), 'LineWidth', 5, 'Color', 'k');
xlim(cueOnsetLfpXLim);
xlabel('Time from Cue Onset (ms)');
title(sprintf('Cue Evoked Potential, Raw and Average (N=%d)', size(lfpAroundCueOnset, 2)));

%%
figure;
hold on;

lfpAroundCueOnsetByLoc = cell(nLoc, 1);
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundCueOnsetByLoc{i} = createdatamatc(adjLfp, cueOnsetByLoc{i}, Fs, cueOnsetLfpWindow);
    plot(cueOnsetLfpT, mean(lfpAroundCueOnsetByLoc{i}, 2), 'LineWidth', 2);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(cueOnsetByLoc{i}));
end

origYLim = ylim();
plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(cueOnsetLfpXLim);
xlabel('Time from Cue Onset (ms)');
title(sprintf('Cue Evoked Potential by Location (N=%d)', size(lfpAroundCueOnset, 2)));
legend(legendEntry, ...
        'Location', 'NorthWest');

%%
arrayOnsetLfpWindow = [1.5 1];
arrayOnsetLfpXLim = [-0.7 0.6] * Fs;
arrayOnsetLfpT = -arrayOnsetLfpWindow(1) * Fs : arrayOnsetLfpWindow(2) * Fs - 1;
lfpAroundArrayOnset = createdatamatc(adjLfp, arrayOnset, Fs, arrayOnsetLfpWindow);
assert(~any(any(isnan(lfpAroundArrayOnset))));

figure;
hold on;
for i = 1:size(lfpAroundArrayOnset, 2)
    plot(arrayOnsetLfpT, lfpAroundArrayOnset(:,i));
end

plot(arrayOnsetLfpT, mean(lfpAroundArrayOnset, 2), 'LineWidth', 5, 'Color', 'k');
xlim(arrayOnsetLfpXLim);
xlabel('Time from Array Onset (ms)');
title(sprintf('Array Evoked Potential (N=%d)', size(lfpAroundArrayOnset, 2)));

%%
figure;
hold on;

arrayOnsetLfpXLim = [-0.7 1] * Fs;
lfpAroundArrayOnsetHoldByLoc = cell(nLoc, 1);
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundArrayOnsetHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
    plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetHoldByLoc{i}, 2), 'LineWidth', 2);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetHoldByLoc{i}));
end

origYLim = ylim();
plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(arrayOnsetLfpXLim);
xlabel('Time from Array Onset (ms)');
title(sprintf('Array Evoked Potential by Location -- Hold Target (N=%d)', numel(arrayOnsetHold)));
legend(legendEntry, ...
        'Location', 'NorthWest');

%%
figure;
hold on;

arrayOnsetLfpXLim = [-0.7 1] * Fs;
lfpAroundArrayOnsetShortHoldByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetShortHoldByLocFilt = cell(nLoc, 1);
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundArrayOnsetShortHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetShortHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
    assert(~any(any(isnan(lfpAroundArrayOnsetShortHoldByLoc{i}))));
    lfpAroundArrayOnsetShortHoldByLocFilt{i} = filtfilt(bFirLowPass, 1, lfpAroundArrayOnsetShortHoldByLoc{i});
    plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetShortHoldByLocFilt{i}, 2), 'LineWidth', 2);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetShortHoldByLoc{i}));
end

origYLim = ylim();
plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(arrayOnsetLfpXLim);
xlabel('Time from Array Onset (ms)');
title(sprintf('Array Evoked Potential by Location -- Hold Target Short Dur (N=%d)', nTrialShortHoldDur));
legend(legendEntry, ...
        'Location', 'NorthWest');

%%
figure;
hold on;

arrayOnsetLfpXLim = [-0.7 1] * Fs;

lfpAroundArrayOnsetLongHoldByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetLongHoldByLocFilt = cell(nLoc, 1);
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundArrayOnsetLongHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetLongHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
    assert(~any(any(isnan(lfpAroundArrayOnsetLongHoldByLoc{i}))));
    lfpAroundArrayOnsetLongHoldByLocFilt{i} = filtfilt(bFirLowPass, 1, lfpAroundArrayOnsetLongHoldByLoc{i});
    plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetLongHoldByLocFilt{i}, 2), 'LineWidth', 2);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetLongHoldByLoc{i}));
end

origYLim = ylim();
plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(arrayOnsetLfpXLim);
xlabel('Time from Array Onset (ms)');
title(sprintf('Array Evoked Potential by Location -- Hold Target Long Dur (N=%d)', nTrialLongHoldDur));
legend(legendEntry, ...
        'Location', 'NorthWest');

%%
figure;
hold on;

arrayOnsetLfpXLim = [-0.7 0.6] * Fs;
lfpAroundArrayOnsetRelByLoc = cell(nLoc, 1);
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundArrayOnsetRelByLoc{i} = createdatamatc(adjLfp, arrayOnsetRelByLoc{i}, Fs, arrayOnsetLfpWindow);
    plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetRelByLoc{i}, 2), 'LineWidth', 2);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetRelByLoc{i}));
end

origYLim = ylim();
plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(arrayOnsetLfpXLim);
xlabel('Time from Array Onset (ms)');
title(sprintf('Array Evoked Potential by Location -- Release Target (N=%d)', numel(arrayOnsetRel)));
legend(legendEntry, ...
        'Location', 'NorthWest');
    
%%
targetDimLfpWindow = [1.5 1];
targetDimLfpXLim = [-0.7 0.4] * Fs;
targetDimLfpT = -targetDimLfpWindow(1) * Fs : targetDimLfpWindow(2) * Fs - 1;
lfpAroundTargetDim = createdatamatc(adjLfp, targetDim, Fs, targetDimLfpWindow);
assert(~any(any(isnan(lfpAroundTargetDim))));

figure;
hold on;
for i = 1:size(lfpAroundTargetDim, 2)
    plot(arrayOnsetLfpT, lfpAroundTargetDim(:,i));
end

plot(targetDimLfpT, mean(lfpAroundTargetDim, 2), 'LineWidth', 5, 'Color', 'k');
xlim(targetDimLfpXLim);
xlabel('Time from Target Dimming (ms)');
title(sprintf('Target Dimming Potential (N=%d)', size(lfpAroundTargetDim, 2)));

%%
figure;
hold on;

targetDimLfpXLim = [-0.7 0.4] * Fs;
lfpAroundTargetDimByLoc = cell(nLoc, 1);
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundTargetDimByLoc{i} = createdatamatc(adjLfp, targetDimByLoc{i}, Fs, targetDimLfpWindow);
    plot(targetDimLfpT, mean(lfpAroundTargetDimByLoc{i}, 2), 'LineWidth', 2);
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetHoldByLoc{i}));
end

origYLim = ylim();
plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(targetDimLfpXLim);
xlabel('Time from Target Dimming (ms)');
title(sprintf('Target Dimming Evoked Potential by Location (N=%d)', numel(targetDim)));
legend(legendEntry, ...
        'Location', 'NorthWest');

%%

load('params.mat');
params.fpass = [10 40];

specGramWindow = [0.2 0.1];
[S,t,f,SErr] = mtspecgramc(lfpAroundArrayOnsetHoldByLoc{1}, specGramWindow, params);

figure;
hold on;
plot_matrix(S, t - arrayOnsetLfpWindow(1), f, 'l');
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(params.fpass);

%%
specGramWindow = [0.2 0.1];
[S,t,f,SErr] = mtspecgramc(lfpAroundCueOnsetByLoc{1}, specGramWindow, params);

figure;
hold on;
plot_matrix(S, t - cueOnsetLfpWindow(1), f, 'l');
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(params.fpass);

%%
params.fpass = [4 40];
arrayOnsetLfpSelection = arrayOnsetLfpT >= 400 & arrayOnsetLfpT <= 1000;
[S,f,SErr] = mtspectrumc(lfpAroundArrayOnsetHoldByLoc{1}(arrayOnsetLfpSelection,:), params);

figure;
hold on;
plot_vector(S, f, 'l');
xlim(params.fpass);

%%
params.fpass = [4 40];
arrayOnsetLfpSelection = arrayOnsetLfpT >= -600 & arrayOnsetLfpT <= 0;
[S,f,SErr] = mtspectrumc(lfpAroundArrayOnsetHoldByLoc{4}(arrayOnsetLfpSelection,:), params);

figure;
hold on;
plot_vector(S, f, 'l');
xlim(params.fpass);

%%