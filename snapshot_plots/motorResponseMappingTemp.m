% function spikeVEPMapping(sessionInd)
% saccade response (immediately pre/post task only)
% x = 952-966, y = 961-976
% lever response around 164.2 +/- 1 (or -2500 or 2500)
% juice response
% for now, use VEPM blocks

% plot saccade and lever response aligned to juice event

clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';

sessionInd = 8;
doPlot = 1;

%% load recording information
recordingInfo = readRecordingInfo();
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% setup and load data
fprintf('\n-------------------------------------------------------\n');
fprintf('Juice and Lever Response Mapping Analysis - Spikes\n');
fprintf('Loading %s...\n', pl2FilePath);

tic;
isLoadSpikes = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
isLoadDirect = 1;
D = loadPL2(pl2FilePath, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad); 

processedDataDir = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

blockName = strjoin(blockNames(gratingsTask3DIndices), '-');

%% remove spike and event times not during RFM task to save memory
D = trimSpikeTimesAndEvents(D, gratingsTask3DIndices);

%% clean up juice events
% store only the first juice event
juiceEvent = D.events{8};
firstJuiceEvent = nan(1, 1);
count = 0;
diffJuiceEvent = diff(juiceEvent);
maxDiffJuiceEventSameTrial = 1; % seconds
for i = 1:numel(juiceEvent)
    if i >= 2 && diffJuiceEvent(i-1) < maxDiffJuiceEventSameTrial
        continue;
    end
    count = count + 1;
    firstJuiceEvent(count) = juiceEvent(i);
end

% for VEPM mapping
if strcmp(sessionName, 'M20170127') || strcmp(sessionName, 'M20170130') || strcmp(sessionName, 'M20170201')
    trialStartEvent = D.events{5};
else
    trialStartEvent = D.events{2};
end

% for gratings task -- EVT04 is cue onset
trialStartEvent = D.events{4};

%%
direct1LockedToJuice = createdatamatc(D.adjDirects(1,:)', firstJuiceEvent, 1000, [1 1]);
direct2LockedToJuice = createdatamatc(D.adjDirects(2,:)', firstJuiceEvent, 1000, [1 1]);
direct3LockedToJuice = createdatamatc(D.adjDirects(3,:)', firstJuiceEvent, 1000, [1 1]);

figure_tr_inch(13, 6);
subaxis(1, 3, 1);
plot(direct1LockedToJuice);
ylim(959 + [-25 25]);
subaxis(1, 3, 2);
plot(direct2LockedToJuice);
ylim(968.5 + [-25 25]);
subaxis(1, 3, 3);
plot(direct3LockedToJuice);
ylim(164.2 + [-3 3]);

direct1LockedToTrialStart = createdatamatc(D.adjDirects(1,:)', trialStartEvent, 1000, [1 1]);
direct2LockedToTrialStart = createdatamatc(D.adjDirects(2,:)', trialStartEvent, 1000, [1 1]);
direct3LockedToTrialStart = createdatamatc(D.adjDirects(3,:)', trialStartEvent, 1000, [1 1]);

figure_tr_inch(13, 6);
subaxis(1, 3, 1);
plot(direct1LockedToTrialStart);
ylim(959 + [-25 25]);
subaxis(1, 3, 2);
plot(direct2LockedToTrialStart);
ylim(968.5 + [-25 25]);
subaxis(1, 3, 3);
plot(direct3LockedToTrialStart);
ylim(164.2 + [-3 3]);

%% 
% need to pass cueLoc and nLoc...
% note that in main code, only correct trialStartEvents are used
fixationAndLeverTimes = getFixationAndLeverTimes(D, trialStartEvent, firstJuiceEvent, 0, 0);
struct2var(fixationAndLeverTimes);

%%
figure;
hold on;
plot(firstLeverPressTimesPreCue, zeros(size(firstLeverPressTimesPreCue)), '.', 'MarkerSize', 20);
plot(firstEnterFixationTimesPreCue, ones(size(firstEnterFixationTimesPreCue)), '.', 'MarkerSize', 20);
plot(trialStartEvent, 2*ones(size(trialStartEvent)), '.', 'MarkerSize', 20);
plot(firstLeverReleaseTimesAroundJuice, 3*ones(size(firstLeverReleaseTimesAroundJuice)), '.', 'MarkerSize', 20);
plot(firstJuiceEvent, 4*ones(size(firstJuiceEvent)), '.', 'MarkerSize', 20);
plot(firstExitFixationTimesAroundJuice, 5*ones(size(firstExitFixationTimesAroundJuice)), '.', 'MarkerSize', 20);
% xlim([0 300]);

% 0 = lever press before trial start
% 1 = enter fixation before trial start (could be before or after 0)
% 2 = trial start
% 3 = lever release before trial start
% 4 = juice event
% 5 = exit fixation before trial start (could be before or after 3 or 4)
% 0 and 1 must precede 2
% 3, 4, 5 may not happen
% 3 must precede 4
% 5 must be around 4

%% check event numbers

% the same fixation or lever press event may be associated with multiple
% trial start events (rare)
% there may be lost events due to edge effects as well
fprintf('\n');
fprintf('%d trial start events.\n', numel(trialStartEvent));
fprintf('%d first lever press times preceding trial start.\n', numel(firstLeverPressTimesPreCue));
fprintf('%d first enter fixation times preceding trial start.\n', numel(firstEnterFixationTimesPreCue));
fprintf('\n');
fprintf('%d first juice events.\n', numel(firstJuiceEvent));
fprintf('%d first lever release times preceding first juice.\n', numel(firstLeverReleaseTimesAroundJuice));
fprintf('%d first exit fixation times around first juice.\n', numel(firstExitFixationTimesAroundJuice));
fprintf('\n');
fprintf('%d other lever press times.\n', numel(otherLeverPressTimes));
fprintf('%d other lever release times.\n', numel(otherLeverReleaseTimes));

stop

%% animate fixation during a trial (use trialStartEvent)
% alternatively, center on juice event firstJuiceEvent
eventTimes = trialStartEvent;

for i = 1:numel(eventTimes)
    figure_tr_inch(10, 10); hold on;
    set(gcf, 'Color', 'w');
    plot(voltages.fixationXVoltage + [-50 50], [voltages.fixationYVoltage voltages.fixationYVoltage], 'Color', 0.3*ones(3, 1));
    plot([voltages.fixationXVoltage voltages.fixationXVoltage], voltages.fixationYVoltage + [-50 50], 'Color', 0.3*ones(3, 1));
    
    % draw a circle at fixation
    diameter = voltages.rangeFixationXVoltage * 2;
    theta = linspace(0, 2 * pi, 100);
    rho = ones(1, 100) * diameter / 2;
    [circX,circY] = pol2cart(theta, rho);
    circX = circX + voltages.fixationXVoltage;
    circY = circY + voltages.fixationYVoltage;
    plot(circX, circY, 'Color', 0.3 * ones(3, 1));
    
    xlim(voltages.fixationXVoltage + [-50 50]);
    ylim(voltages.fixationYVoltage + [-50 50]);
    axis square;
    
    for j = -1000:10:3000
        tInd = round(eventTimes(i) * D.directFs + j);
        plot(D.adjDirects(1,tInd), D.adjDirects(2,tInd), ...
                '.', 'MarkerSize', 20, 'Color', rand(3, 1));
        if abs(D.adjDirects(3,tInd) - voltages.leverPressedVoltage) > voltages.rangeLeverPressedVoltage
            set(gca, 'Color', [1 0.5 0.5]);
        else
            set(gca, 'Color', [0.5 1 0.5]);
        end
        title(sprintf('Trial %d: t = %d ms', i, j));
        pause(0.05);
    end
    pause;
    close;
end

%% show event locked psths for enter/exit fixation and press/rel lever
eventAlignmentWindow = [0.6 0.6]; % secs before, secs after
psthWindowOffset = [-0.5 0.5]; % secs before, secs after
kernelSigma = 0.01;
psthT = computeTForSpdf(eventAlignmentWindow(1), psthWindowOffset, kernelSigma);
t = psthT - eventAlignmentWindow(1);

minAbsNormMaxRespJuiceDriven = 3; % normed value is z-scored (SD units)
minProportionJuiceEventsWithSpikeJuiceDriven = 0.1;

leverArtifactTestWindow = [-0.01 0.01];
leverArtifactTLogical = getTimeLogicalWithTolerance(t, leverArtifactTestWindow);

leverArtifactBaselineWindow = [-0.325 -0.025];
leverArtifactBaselineLogical = getTimeLogicalWithTolerance(t, leverArtifactBaselineWindow);
nSDOverTimeOverMeanBaselineForArtifact = 3;

nUnits = numel(D.allSpikeStructs);
for i = 1:nUnits
    alignedSpikeTimesEnter = createdatamatpt(D.allSpikeStructs{i}.ts, firstEnterFixationTimesPreCue, eventAlignmentWindow);
    psthResponseEnter = fixedPsth(alignedSpikeTimesEnter, kernelSigma, 2, psthT);

    alignedSpikeTimesExit = createdatamatpt(D.allSpikeStructs{i}.ts, firstExitFixationTimesAroundJuice, eventAlignmentWindow);
    psthResponseExit = fixedPsth(alignedSpikeTimesExit, kernelSigma, 2, psthT);
    
    alignedSpikeTimesOtherEnter = createdatamatpt(D.allSpikeStructs{i}.ts, otherEnterFixationTimes, eventAlignmentWindow);
    psthResponseOtherEnter = fixedPsth(alignedSpikeTimesOtherEnter, kernelSigma, 2, psthT);

    alignedSpikeTimesOtherExit = createdatamatpt(D.allSpikeStructs{i}.ts, otherExitFixationTimes, eventAlignmentWindow);
    psthResponseOtherExit = fixedPsth(alignedSpikeTimesOtherExit, kernelSigma, 2, psthT);
    
    alignedSpikeTimesPress = createdatamatpt(D.allSpikeStructs{i}.ts, firstLeverPressTimesPreCue, eventAlignmentWindow);
    psthResponsePress = fixedPsth(alignedSpikeTimesPress, kernelSigma, 2, psthT);

    alignedSpikeTimesRelease = createdatamatpt(D.allSpikeStructs{i}.ts, firstLeverReleaseTimesAroundJuice, eventAlignmentWindow);
    psthResponseRelease = fixedPsth(alignedSpikeTimesRelease, kernelSigma, 2, psthT);
    
    alignedSpikeTimesOtherPress = createdatamatpt(D.allSpikeStructs{i}.ts, otherLeverPressTimes, eventAlignmentWindow);
    psthResponseOtherPress = fixedPsth(alignedSpikeTimesOtherPress, kernelSigma, 2, psthT);

    alignedSpikeTimesOtherRelease = createdatamatpt(D.allSpikeStructs{i}.ts, otherLeverReleaseTimes, eventAlignmentWindow);
    psthResponseOtherRelease = fixedPsth(alignedSpikeTimesOtherRelease, kernelSigma, 2, psthT);
    
    % if there is a peak in spiking in the +/- 10 ms window around other 
    % press or other release compared to the 300 ms before, then
    % the unit contains lever artifacts OR the unit is responsive to lever
    % presses/releases
    peakAroundOtherPress = max(psthResponseOtherPress(leverArtifactTLogical));
    peakAroundOtherRelease = max(psthResponseOtherRelease(leverArtifactTLogical));
    
    meanLeverArtifactBaselineOtherPress = mean(psthResponseOtherPress(leverArtifactBaselineLogical));
    sdOverTimeLeverArtifactBaselineOtherPress = std(psthResponseOtherPress(leverArtifactBaselineLogical));
    leverArtifactThreshOtherPress = meanLeverArtifactBaselineOtherPress + ...
            nSDOverTimeOverMeanBaselineForArtifact * sdOverTimeLeverArtifactBaselineOtherPress;
        
    meanLeverArtifactBaselineOtherRelease = mean(psthResponseOtherRelease(leverArtifactBaselineLogical));
    sdOverTimeLeverArtifactBaselineOtherRelease = std(psthResponseOtherRelease(leverArtifactBaselineLogical));
    leverArtifactThreshOtherRelease = meanLeverArtifactBaselineOtherRelease + ...
            nSDOverTimeOverMeanBaselineForArtifact * sdOverTimeLeverArtifactBaselineOtherRelease;
        
    if peakAroundOtherPress > leverArtifactThreshOtherPress || ...
            peakAroundOtherRelease > leverArtifactThreshOtherRelease
        hasLeverArtifact = true;
    else
        hasLeverArtifact = false;
    end
    
    cols = lines(4);
    figure_tr_inch(10, 10);
    hold on;
    plot(t, psthResponseEnter, '-', 'LineWidth', 2, 'Color', cols(1,:)); % blue
    plot(t, psthResponseExit, '-', 'LineWidth', 2, 'Color', cols(2,:)); % red
    plot(t, psthResponsePress, '-', 'LineWidth', 2, 'Color', cols(3,:)); % yellow
    plot(t, psthResponseRelease, '-', 'LineWidth', 2, 'Color', cols(4,:)); % purple
    plot(t, psthResponseOtherEnter, ':', 'LineWidth', 2, 'Color', cols(1,:)); % blue
    plot(t, psthResponseOtherExit, ':', 'LineWidth', 2, 'Color', cols(2,:)); % red
    plot(t, psthResponseOtherPress, ':', 'LineWidth', 2, 'Color', cols(3,:)); % yellow
    plot(t, psthResponseOtherRelease, ':', 'LineWidth', 2, 'Color', cols(4,:)); % purple
    
    origYLim = ylim();
    plot([0 0], ylim, '-', 'Color', 0.3*ones(3, 1));
    ylim(origYLim);
    
    legend({'Enter Fixation', 'Exit Fixation', 'Press Lever', 'Release Lever'}, ...
            'Location', 'NorthWest');
    
    if hasLeverArtifact
        set(gca, 'Color', [1 0.5 0.5]);
    else
        set(gca, 'Color', [0.5 1 0.5]);
    end
    set(gca, 'FontSize', 14);
    
    title(sprintf('Unit %d', i));
    drawnow;
%     pause;
%     close;
end

% note: post exit saccade and release lever, there can be a major juice /
% licking artifact, though technically he could be licking just before too