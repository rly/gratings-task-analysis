function lfpAnalysis(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, lfpChannelsToLoad, isZeroDistractors)
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

v = 11;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task LFP Analysis\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('LFP Channel to Load: %d\n', lfpChannelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Is zero distractors: %s\n', isZeroDistractors);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

nLoc = 4;

%% input check
% assert(numel(lfpChannelsToLoad) == 1);
assert(numel(lfpChannelsToLoad) > 1);
assert(numel(lfpChannelsToLoad) <= 32);

%% load recording information

if isZeroDistractors
    scriptName = 'LFP_GRATINGS_0D';
    taskName = 'GRATINGS_0D';
else
    scriptName = 'LFP_GRATINGS';
    taskName = 'GRATINGS';
end

[R, D, processedDataDir, blockName] = loadRecordingData(processedDataRootDir, ...
        dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, lfpChannelsToLoad, ...
        taskName, scriptName, 1, 1);
sessionName = R.sessionName;
areaName = R.areaName;

fprintf('Processing %s...\n', sessionName);

dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
gratingsTaskLogDir = sprintf('%s/%s', dataDir, sessionName(2:end));

% process events and sort them into different conditions
UE = getUsefulEvents2(gratingsTaskLogDir, R.gratingsTaskLogIndices, 4, D, blockName);

fileNamePrefix = sprintf('%s-ch%d-ch%d-%s', sessionName, lfpChannelsToLoad([1 end]), blockName);

%% preprocess LFPs
tic;
Fs = D.lfpFs;
nChannels = D.nLfpCh;

D.adjLfpsClean = interpolateLfpOverSpikeTimes(D.adjLfps, lfpChannelsToLoad, Fs, D.allMUAStructs);

hiCutoffFreq = 100;
[channelDataNorm,isTrialOutlier,isNoisyChannel,meanBeforeCAR] = preprocessLfpsGratingsTask(D.adjLfpsClean, Fs, D.lfpNames, UE, ...
        processedDataDir, fileNamePrefix, hiCutoffFreq, v); % TODO USE isTrialOutlier
D.adjLfps = [];
D.adjLfpsClean = [];
D.adjDirects = [];

% LAST CHANNEL IS MEAN BEFORE COMMON AVERAGE REFERENCING
channelDataNorm(size(channelDataNorm, 1)+1,:) = meanBeforeCAR;

fprintf('Took %0.1f minutes.\n', toc/60);

%% save evoked file
fprintf('Extracting event-locked LFPs...\n');
tic;
saveFileName = sprintf('%s/%s-evokedLfps-v%d.mat', ...
        processedDataDir, fileNamePrefix, v);
computeEvokedLfps(saveFileName, channelDataNorm, Fs, nLoc, UE);
fprintf('Took %0.1f minutes.\n', toc/60);

%% load evoked file
% clear channelDataNorm;
saveFileName = sprintf('%s/%s-evokedLfps-v%d.mat', ...
        processedDataDir, fileNamePrefix, v);
EL = load(saveFileName);

%%
% save('temp_workspace.mat');

%%
% clear;
% load('temp_workspace.mat');

%% combined line plot
cols = lines(nLoc);
yScale = 1;
ySep = -1;
xBounds = [-0.4 0.4];
yBounds = [(nChannels+4)*ySep -1*ySep];

baselineWindowOffset = [-0.3 0];
baselineInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, baselineWindowOffset);
preCueBaseline = nan(nChannels+1, nLoc);
meanInd = nChannels+1;
for j = 1:nChannels+1    
    for m = 1:nLoc
        preCueBaseline(j,m) = squeeze(mean(mean(EL.cueOnsetLfp.lfp(j,UE.cueLoc == m,baselineInd), 2), 3));
    end
end

figure_tr_inch(20, 10);
subaxis(1, 5, 1, 'SH', 0.02, 'MT', 0.06);
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
for j = 1:nChannels
    for m = 1:nLoc          
        plot(EL.enterFixationLfp.t, (squeeze(mean(EL.enterFixationLfp.lfp(j,UE.cueLoc == m,:), 2)) - preCueBaseline(j,m))*yScale + j*ySep, 'Color', cols(m,:), 'LineWidth', 1);
    end
end
for m = 1:nLoc          
    plot(EL.enterFixationLfp.t, (squeeze(mean(EL.enterFixationLfp.lfp(meanInd,UE.cueLoc == m,:), 2)) - preCueBaseline(meanInd,m))*yScale + (meanInd+1)*ySep, 'Color', cols(m,:), 'LineWidth', 2);
end
xlim(xBounds);
ylim(yBounds);
set(gca, 'YTick', [-34 -32:-1]);
set(gca, 'YTickLabel', ['CAR' strsplit(num2str(flip(1:32)))]);
set(gca, 'box', 'off');
xlabel('Time from Start Fixation (s)');
ylabel('Channel Number');
title('Start Fixation');

subaxis(1, 5, 2);
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
for j = 1:nChannels    
    for m = 1:nLoc
        plot(EL.cueOnsetLfp.t, (squeeze(mean(EL.cueOnsetLfp.lfp(j,UE.cueLoc == m,:), 2)) - preCueBaseline(j,m))*yScale + j*ySep, 'Color', cols(m,:), 'LineWidth', 1);
    end
end
for m = 1:nLoc
    plot(EL.cueOnsetLfp.t, (squeeze(mean(EL.cueOnsetLfp.lfp(meanInd,UE.cueLoc == m,:), 2)) - preCueBaseline(meanInd,m))*yScale + (meanInd+1)*ySep, 'Color', cols(m,:), 'LineWidth', 2);
end
xlim(xBounds);
ylim(yBounds);
set(gca, 'YTick', [-34 -32:-1]);
set(gca, 'YTickLabel', ['CAR' strsplit(num2str(flip(1:32)))]);
set(gca, 'box', 'off');
xlabel('Time from Cue Onset (s)');
title('Cue Onset');

subaxis(1, 5, 3);
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
for j = 1:nChannels    
    for m = 1:nLoc
        plot(EL.arrayOnsetHoldBalLfp.t, (squeeze(mean(EL.arrayOnsetHoldBalLfp.lfp(j,UE.cueLocHoldBal == m,:), 2)) - preCueBaseline(j,m))*yScale + j*ySep, 'Color', cols(m,:), 'LineWidth', 1);
    end
end
for m = 1:nLoc
    plot(EL.arrayOnsetHoldBalLfp.t, (squeeze(mean(EL.arrayOnsetHoldBalLfp.lfp(meanInd,UE.cueLocHoldBal == m,:), 2)) - preCueBaseline(meanInd,m))*yScale + (meanInd+1)*ySep, 'Color', cols(m,:), 'LineWidth', 2);
end
xlim(xBounds);
ylim(yBounds);
set(gca, 'YTick', [-34 -32:-1]);
set(gca, 'YTickLabel', ['CAR' strsplit(num2str(flip(1:32)))]);
set(gca, 'box', 'off');
xlabel('Time from Array Onset (s)');
title('Array Onset');

subaxis(1, 5, 4);
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
for j = 1:nChannels    
    for m = 1:nLoc
        plot(EL.targetDimBalLfp.t, (squeeze(mean(EL.targetDimBalLfp.lfp(j,UE.cueLocHoldBal == m,:), 2)) - preCueBaseline(j,m))*yScale + j*ySep, 'Color', cols(m,:), 'LineWidth', 1);
    end
end
for m = 1:nLoc
    plot(EL.targetDimBalLfp.t, (squeeze(mean(EL.targetDimBalLfp.lfp(meanInd,UE.cueLocHoldBal == m,:), 2)) - preCueBaseline(meanInd,m))*yScale + (meanInd+1)*ySep, 'Color', cols(m,:), 'LineWidth', 2);
end
xlim(xBounds);
ylim(yBounds);
set(gca, 'YTick', [-34 -32:-1]);
set(gca, 'YTickLabel', ['CAR' strsplit(num2str(flip(1:32)))]);
set(gca, 'box', 'off');
xlabel('Time from Target Dimming (s)');
title('Target Dimming');

subaxis(1, 5, 5);
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
for j = 1:nChannels    
    for m = 1:nLoc
        plot(EL.exitFixationLfp.t, (squeeze(mean(EL.exitFixationLfp.lfp(j,UE.cueLoc == m,:), 2)) - preCueBaseline(j,m))*yScale + j*ySep, 'Color', cols(m,:), 'LineWidth', 1);
    end
end
for m = 1:nLoc
    plot(EL.exitFixationLfp.t, (squeeze(mean(EL.exitFixationLfp.lfp(meanInd,UE.cueLoc == m,:), 2)) - preCueBaseline(meanInd,m))*yScale + (meanInd+1)*ySep, 'Color', cols(m,:), 'LineWidth', 2);
end
xlim(xBounds);
ylim(yBounds);
set(gca, 'YTick', [-34 -32:-1]);
set(gca, 'YTickLabel', ['CAR' strsplit(num2str(flip(1:32)))]);
set(gca, 'box', 'off');
xlabel('Time from Exit Fixation (s)');
title('Exit Fixation');

axBig = axes('Position', [0.04 0.035 0.92 0.93], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on');
title(axBig, sprintf('Session %s (Ch %d-%d)', sessionName, lfpChannelsToLoad([1 end])), 'FontSize', 14);

plotFileName = sprintf('%s/%s-allFP-evokedLfps-combLinePlot-v%d.png', ...
        processedDataDir, fileNamePrefix, v);
fprintf('Saving figure to file %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% separated line plot
cols = lines(4);
yScale = 1;
ySep = -1;
xBounds = [-0.4 0.4];
yBounds = [(nChannels+2)*ySep -1*ySep];

baselineWindowOffset = [-0.3 0];
baselineInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, baselineWindowOffset);

figure_tr_inch(20, 10);
clf;
set(gcf, 'renderer', 'painters');

for m = 1:nLoc
    
preCueBaseline = squeeze(mean(mean(EL.cueOnsetLfp.lfp(:,UE.cueLoc == m,baselineInd), 2), 3));

subaxis(nLoc, 5, (m-1)*5 + 1, 'SH', 0.02, 'MT', 0.06);
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
for j = 1:nChannels    
        plot(EL.enterFixationLfp.t, (squeeze(mean(EL.enterFixationLfp.lfp(j,UE.cueLoc == m,:), 2)) - preCueBaseline(j))*yScale + j*ySep, 'Color', cols(m,:));
end
xlim(xBounds);
ylim(yBounds);
set(gca, 'box', 'off');
if m == 1
    title('Start Fixation');
end
if m == nLoc
    xlabel('Time from Start Fixation (s)');
end
ylabel('Channel Number');
set(gca, 'YTick', -30:5:0);
set(gca, 'YTickLabel', strsplit(num2str(flip(0:5:30))));

subaxis(nLoc, 5, (m-1)*5 + 2);
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
for j = 1:nChannels    
        plot(EL.cueOnsetLfp.t, (squeeze(mean(EL.cueOnsetLfp.lfp(j,UE.cueLoc == m,:), 2)) - preCueBaseline(j))*yScale + j*ySep, 'Color', cols(m,:));
end
xlim(xBounds);
ylim(yBounds);
set(gca, 'box', 'off');
if m == 1
    title('Cue Onset');
end
if m == nLoc
    xlabel('Time from Cue Onset (s)');
end
set(gca, 'YTick', -30:5:0);
set(gca, 'YTickLabel', strsplit(num2str(flip(0:5:30))));

subaxis(nLoc, 5, (m-1)*5 + 3);
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
for j = 1:nChannels    
        plot(EL.arrayOnsetHoldBalLfp.t, (squeeze(mean(EL.arrayOnsetHoldBalLfp.lfp(j,UE.cueLocHoldBal == m,:), 2)) - preCueBaseline(j))*yScale + j*ySep, 'Color', cols(m,:));
end
xlim(xBounds);
ylim(yBounds);
set(gca, 'box', 'off');
if m == 1
    title('Array Onset');
end
if m == nLoc
    xlabel('Time from Array Onset (s)');
end
set(gca, 'YTick', -30:5:0);
set(gca, 'YTickLabel', strsplit(num2str(flip(0:5:30))));

subaxis(nLoc, 5, (m-1)*5 + 4);
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
for j = 1:nChannels    
        plot(EL.targetDimBalLfp.t, (squeeze(mean(EL.targetDimBalLfp.lfp(j,UE.cueLocHoldBal == m,:), 2)) - preCueBaseline(j))*yScale + j*ySep, 'Color', cols(m,:));
end
xlim(xBounds);
ylim(yBounds);
set(gca, 'box', 'off');
if m == 1
    title('Target Dimming');
end
if m == nLoc
    xlabel('Time from Target Dimming (s)');
end
set(gca, 'YTick', -30:5:0);
set(gca, 'YTickLabel', strsplit(num2str(flip(0:5:30))));

subaxis(nLoc, 5, (m-1)*5 + 5);
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
for j = 1:nChannels    
        plot(EL.exitFixationLfp.t, (squeeze(mean(EL.exitFixationLfp.lfp(j,UE.cueLoc == m,:), 2)) - preCueBaseline(j))*yScale + j*ySep, 'Color', cols(m,:));
end
xlim(xBounds);
ylim(yBounds);
set(gca, 'box', 'off');
if m == 1
    title('Exit Fixation');
end
if m == nLoc
    xlabel('Time from Exit Fixation (s)');
end
set(gca, 'YTick', -30:5:0);
set(gca, 'YTickLabel', strsplit(num2str(flip(0:5:30))));

end % end for each location

axBig = axes('Position', [0.04 0.035 0.92 0.93], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on');
title(axBig, sprintf('Session %s (Ch %d-%d)', sessionName, lfpChannelsToLoad([1 end])), 'FontSize', 14);

plotFileName = sprintf('%s/%s-allFP-evokedLfps-sepLinePlot-v%d.png', ...
        processedDataDir, fileNamePrefix, v);
fprintf('Saving figure to file %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% image plot
xBounds = [-0.4 0.4];
yBounds = [1 nChannels];
cBounds = [-1 1];

baselineWindowOffset = [-0.3 0];
baselineInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, baselineWindowOffset);

figure_tr_inch(20, 10); 
clf;
set(gcf, 'renderer', 'painters');

for m = 1:nLoc
    preCueBaseline = squeeze(mean(mean(EL.cueOnsetLfp.lfp(:,UE.cueLoc == m,baselineInd), 2), 3));

    subaxis(nLoc, 5, (m-1)*5 + 1, 'SH', 0.02, 'MT', 0.06);
    hold on;
    imagesc(EL.enterFixationLfp.t, 1:nChannels, squeeze(mean(EL.enterFixationLfp.lfp(:,UE.cueLoc == m,:), 2)) - preCueBaseline);
    plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
    set(gca, 'YDir', 'reverse');
    xlim(xBounds);
    ylim(yBounds);
    caxis(cBounds);
    set(gca, 'box', 'off');
    if m == 1
        title('Start Fixation');
    end
    if m == nLoc
        xlabel('Time from Start Fixation (s)');
    end
    ylabel('Channel Number');

    subaxis(nLoc, 5, (m-1)*5 + 2);
    hold on;
    imagesc(EL.cueOnsetLfp.t, 1:nChannels, squeeze(mean(EL.cueOnsetLfp.lfp(:,UE.cueLoc == m,:), 2)) - preCueBaseline);
    plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
    set(gca, 'YDir', 'reverse');
    xlim(xBounds);
    ylim(yBounds);
    caxis(cBounds);
    set(gca, 'box', 'off');
    if m == 1
        title('Cue Onset');
    end
    if m == nLoc
        xlabel('Time from Cue Onset (s)');
    end

    subaxis(nLoc, 5, (m-1)*5 + 3);
    hold on;
    imagesc(EL.arrayOnsetHoldBalLfp.t, 1:nChannels, squeeze(mean(EL.arrayOnsetHoldBalLfp.lfp(:,UE.cueLocHoldBal == m,:), 2)) - preCueBaseline);
    plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
    set(gca, 'YDir', 'reverse');
    xlim(xBounds);
    ylim(yBounds);
    caxis(cBounds);
    set(gca, 'box', 'off');
    if m == 1
        title('Array Onset');
    end
    if m == nLoc
        xlabel('Time from Array Onset (s)');
    end

    subaxis(nLoc, 5, (m-1)*5 + 4);
    hold on;
    imagesc(EL.targetDimBalLfp.t, 1:nChannels, squeeze(mean(EL.targetDimBalLfp.lfp(:,UE.cueLocHoldBal == m,:), 2)) - preCueBaseline);
    plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
    set(gca, 'YDir', 'reverse');
    xlim(xBounds);
    caxis(cBounds);
    ylim(yBounds);
    set(gca, 'box', 'off');
    if m == 1
        title('Target Dimming');
    end
    if m == nLoc
        xlabel('Time from Target Dimming (s)');
    end

    subaxis(nLoc, 5, (m-1)*5 + 5);
    hold on;
    imagesc(EL.exitFixationLfp.t, 1:nChannels, squeeze(mean(EL.exitFixationLfp.lfp(:,UE.cueLoc == m,:), 2)) - preCueBaseline);
    plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
    set(gca, 'YDir', 'reverse');
    xlim(xBounds);
    ylim(yBounds);
    caxis(cBounds);
    set(gca, 'box', 'off');
    if m == 1
        title('Exit Fixation');
    end
    if m == nLoc
        xlabel('Time from Exit Fixation (s)');
    end
end % end for each location

axBig = axes('Position', [0.04 0.035 0.92 0.93], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on');
title(axBig, sprintf('Session %s (Ch %d-%d)', sessionName, lfpChannelsToLoad([1 end])), 'FontSize', 14);

plotFileName = sprintf('%s/%s-allFP-evokedLfps-sepColorPlot-v%d.png', ...
        processedDataDir, fileNamePrefix, v);
fprintf('Saving figure to file %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');




%% power all channels
cols = lines(6);

baselineWindowOffset = [-0.3 0];
baselineInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, baselineWindowOffset);

cueResponseOffset = [0 0.3];
cueResponseInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, cueResponseOffset);

cueTargetDelayOffset = [-0.3 0];
cueTargetDelayInd = getTimeLogicalWithTolerance(EL.arrayOnsetHoldBalLfp.t, cueTargetDelayOffset);

arrayHoldResponseOffset = [0 0.3];
arrayHoldResponseInd = getTimeLogicalWithTolerance(EL.arrayOnsetHoldBalLfp.t, arrayHoldResponseOffset);

targetDimDelayOffset = [0 0.3];
targetDimDelayInd = getTimeLogicalWithTolerance(EL.targetDimBalLfp.t, targetDimDelayOffset);

params.tapers = [2 3];
params.fpass = [5 60];
params.pad = 2;
params.Fs = D.lfpFs;
params.trialave = 1;
xBounds = params.fpass;
yBounds = [-45 -15];

for channelInd = 1:nChannels+1
    preCueBaselineResponses = cell(nLoc, 1);
    cueResponses = cell(nLoc, 1);
    cueTargetDelayResponses = cell(nLoc, 1);
    arrayOnsetHoldResponses = cell(nLoc, 1);
    targetDimDelayResponses = cell(nLoc, 1);
    
    for m = 1:nLoc
        cueOnsetLfpCurrent = squeeze(EL.cueOnsetLfp.lfp(channelInd,UE.cueLoc == m,:))'; % each column is a trial
        if isempty(cueOnsetLfpCurrent)
            continue;
        end

        arrayOnsetLfpCurrent = squeeze(EL.arrayOnsetLfp.lfp(channelInd,UE.cueLoc == m,:))'; % each column is a trial
        arrayOnsetHoldLfpCurrent = squeeze(EL.arrayOnsetHoldBalLfp.lfp(channelInd,UE.cueLocHoldBal == m,:))'; % each column is a trial
        targetDimLfpCurrent = squeeze(EL.targetDimBalLfp.lfp(channelInd,UE.cueLocHoldBal == m,:))'; % each column is a trial
        
        preCueBaselineResponses{m} = cueOnsetLfpCurrent(baselineInd,:);
        cueResponses{m} = cueOnsetLfpCurrent(cueResponseInd,:);
        cueTargetDelayResponses{m} = arrayOnsetLfpCurrent(cueTargetDelayInd,:);
        arrayOnsetHoldResponses{m} = arrayOnsetHoldLfpCurrent(arrayHoldResponseInd,:);
        targetDimDelayResponses{m} = targetDimLfpCurrent(targetDimDelayInd,:);
    end

    figure_tr_inch(20, 5);
    subaxis(1, 5, 1, 'SH', 0.03);
    hold on;
    for m = 1:nLoc
        if ~isempty(preCueBaselineResponses{m})
            [S,f] = mtspectrumc(preCueBaselineResponses{m}, params);
            plot(f, 10*log10(S), 'LineWidth', 2, 'Color', cols(m,:));
        end
    end
    xlim(xBounds);
    ylim(yBounds);
    title('Pre Cue Baseline');

    subaxis(1, 5, 2);
    hold on;
    for m = 1:nLoc
        if ~isempty(cueResponses{m})
            [S,f] = mtspectrumc(cueResponses{m}, params);
            plot(f, 10*log10(S), 'LineWidth', 2, 'Color', cols(m,:));
        end
    end
    xlim(xBounds);
    ylim(yBounds);
    title('Cue Response');

    subaxis(1, 5, 3);
    hold on;
    for m = 1:nLoc
        if ~isempty(cueTargetDelayResponses{m})
            [S,f] = mtspectrumc(cueTargetDelayResponses{m}, params);
            plot(f, 10*log10(S), 'LineWidth', 2, 'Color', cols(m,:));
        end
    end
    xlim(xBounds);
    ylim(yBounds);
    title('Cue-Target Delay');
    
    subaxis(1, 5, 4);
    hold on;
    for m = 1:nLoc
        if ~isempty(arrayOnsetHoldResponses{m})
            [S,f] = mtspectrumc(arrayOnsetHoldResponses{m}, params);
            plot(f, 10*log10(S), 'LineWidth', 2, 'Color', cols(m,:));
        end
    end
    xlim(xBounds);
    ylim(yBounds);
    title('Array Hold Response');
    
    subaxis(1, 5, 5);
    hold on;
    for m = 1:nLoc
        if ~isempty(targetDimDelayResponses{m})
            [S,f] = mtspectrumc(targetDimDelayResponses{m}, params);
            plot(f, 10*log10(S), 'LineWidth', 2, 'Color', cols(m,:));
        end
    end
    xlim(xBounds);
    ylim(yBounds);
    title('Target-Dim Delay');
    
    suptitle(sprintf('Channel %d', channelInd));
    
    drawnow;
    
    plotFileName = sprintf('%s/%s-FP%03d-power-v%d.png', ...
            processedDataDir, fileNamePrefix, channelInd, v);
    if channelInd == nChannels + 1
        plotFileName = sprintf('%s/%s-FPCAR-power-v%d.png', ...
                processedDataDir, fileNamePrefix, v);
    end
    fprintf('Saving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end

%% look specifically at session M20170130, channel 20 
% right in the middle of the section of putative PI channels with strong
% flash- and cue-evoked potentials

    
%% spectrogram    
periCueOnsetWindowOffset = [-0.4 0.4];
periCueOnsetInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, periCueOnsetWindowOffset);

periArrayOnsetWindowOffset = [-0.4 0.4];
periArrayOnsetInd = getTimeLogicalWithTolerance(EL.arrayOnsetHoldBalLfp.t, periArrayOnsetWindowOffset);

periTargetDimWindowOffset = [-0.4 0.4];
periTargetDimInd = getTimeLogicalWithTolerance(EL.targetDimBalLfp.t, periTargetDimWindowOffset);

params.tapers = [2 3];
params.fpass = [5 60];
params.pad = 2;
params.Fs = D.lfpFs;
params.trialave = 1;
yBounds = params.fpass;
cBounds = [-45 -15];
cDiffBounds = [-1.5 1.5];
movingWin = [0.2 0.05];

for channelInd = 1:nChannels+1
    cueResponses = cell(nLoc, 1);
    arrayOnsetHoldResponses = cell(nLoc, 1);
    targetDimResponses = cell(nLoc, 1);
    
    for m = 1:nLoc
        cueOnsetLfpCurrent = squeeze(EL.cueOnsetLfp.lfp(channelInd,UE.cueLoc == m,:))'; % each column is a trial
        if isempty(cueOnsetLfpCurrent)
            continue;
        end

        arrayOnsetHoldLfpCurrent = squeeze(EL.arrayOnsetHoldBalLfp.lfp(channelInd,UE.cueLocHoldBal == m,:))'; % each column is a trial
        targetDimLfpCurrent = squeeze(EL.targetDimBalLfp.lfp(channelInd,UE.cueLocHoldBal == m,:))'; % each column is a trial
        
        cueResponses{m} = cueOnsetLfpCurrent(periCueOnsetInd,:);
        arrayOnsetHoldResponses{m} = arrayOnsetHoldLfpCurrent(periArrayOnsetInd,:);
        targetDimResponses{m} = targetDimLfpCurrent(periTargetDimInd,:);
    end
    
%     figure_tr_inch(20, 10);
%     
%     for m = 1:nLoc
%         subaxis(nLoc, 3, (m-1)*3+1, 'SH', 0.03);
%         hold on;
%         if ~isempty(cueResponses{m})
%             [S,t,f] = mtspecgramc(cueResponses{m}, movingWin, params);
%             imagesc(t + periCueOnsetWindowOffset(1), f, 10*log10(S'));
%             set(gca, 'YDir', 'normal');
%             xlim(periCueOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
%             ylim(yBounds);
%             caxis(cBounds);
%             xlabel('Time from Cue Onset (s)');
%         end
%     
%         subaxis(nLoc, 3, (m-1)*3+2);
%         hold on;
%         if ~isempty(arrayOnsetHoldResponses{m})
%             [S,t,f] = mtspecgramc(arrayOnsetHoldResponses{m}, movingWin, params);
%             imagesc(t + periArrayOnsetWindowOffset(1), f, 10*log10(S'));
%             set(gca, 'YDir', 'normal');
%             xlim(periArrayOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
%             ylim(yBounds);
%             caxis(cBounds);
%             xlabel('Time from Array Onset (s)');
%         end
%     
%         subaxis(nLoc, 3, (m-1)*3+3);
%         hold on;
%         if ~isempty(targetDimResponses{m})
%             [S,t,f] = mtspecgramc(targetDimResponses{m}, movingWin, params);
%             imagesc(t + periTargetDimWindowOffset(1), f, 10*log10(S'));
%             set(gca, 'YDir', 'normal');
%             xlim(periTargetDimWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
%             ylim(yBounds);
%             caxis(cBounds);
%             xlabel('Time from Target Dimming (s)');
%         end
%     end

%     suptitle(sprintf('Channel %d', channelInd));
    
%     plotFileName = sprintf('%s/%s-FP%03d-powerTfr-v%d.png', ...
%             processedDataDir, fileNamePrefix, channelInd, v);
%     if channelInd == nChannels + 1
%         plotFileName = sprintf('%s/%s-FPCAR-powerTfr-v%d.png', ...
%                 processedDataDir, fileNamePrefix, v);
%     end
%     fprintf('Saving figure to file %s...\n', plotFileName);
%     export_fig(plotFileName, '-nocrop');

    figure_tr_inch(20, 5);
    inRFLoc = 3;
    exRFLoc = 1;
    
    subaxis(1, 3, 1, 'SH', 0.03);
    hold on;
    if ~isempty(cueResponses{inRFLoc})
        [SInRF,t,f] = mtspecgramc(cueResponses{inRFLoc}, movingWin, params);
        [SExRF,t,f] = mtspecgramc(cueResponses{exRFLoc}, movingWin, params);
        SDiff = 10*log10(SInRF') - 10*log10(SExRF');
        imagesc(t + periCueOnsetWindowOffset(1), f, SDiff);
        set(gca, 'YDir', 'normal');
        xlim(periCueOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
        ylim(yBounds);
        caxis(cDiffBounds);
        xlabel('Time from Cue Onset (s)');
        colormap(getCoolWarmMap());
    end

    subaxis(1, 3, 2);
    hold on;
    if ~isempty(arrayOnsetHoldResponses{inRFLoc})
        [SInRF,t,f] = mtspecgramc(arrayOnsetHoldResponses{inRFLoc}, movingWin, params);
        [SExRF,t,f] = mtspecgramc(arrayOnsetHoldResponses{exRFLoc}, movingWin, params);
        SDiff = 10*log10(SInRF') - 10*log10(SExRF');
        imagesc(t + periArrayOnsetWindowOffset(1), f, SDiff);
        set(gca, 'YDir', 'normal');
        xlim(periArrayOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
        ylim(yBounds);
        caxis(cDiffBounds);
        xlabel('Time from Array Onset (s)');
        colormap(getCoolWarmMap());
    end

    subaxis(1, 3, 3);
    hold on;
    if ~isempty(targetDimResponses{inRFLoc})
        [SInRF,t,f] = mtspecgramc(targetDimResponses{inRFLoc}, movingWin, params);
        [SExRF,t,f] = mtspecgramc(targetDimResponses{exRFLoc}, movingWin, params);
        SDiff = 10*log10(SInRF') - 10*log10(SExRF');
        imagesc(t + periTargetDimWindowOffset(1), f, SDiff);
        set(gca, 'YDir', 'normal');
        xlim(periTargetDimWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
        ylim(yBounds);
        caxis(cDiffBounds);
        xlabel('Time from Target Dimming (s)');
        colormap(getCoolWarmMap());
    end
    
    suptitle(sprintf('Channel %d', channelInd));
    
    drawnow;
    plotFileName = sprintf('%s/%s-FP%03d-powerTfrDiff-v%d.png', ...
            processedDataDir, fileNamePrefix, channelInd, v);
    if channelInd == nChannels + 1
        plotFileName = sprintf('%s/%s-FPCAR-powerTfrDiff-v%d.png', ...
                processedDataDir, fileNamePrefix, v);
    end
    fprintf('Saving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end

%% spike field coherence -- chronux first
% units 20a, 20b (45, 46)
% units 19a, 19b, 19c (42, 43, 44)

periCueOnsetWindowOffset = [-0.4 0.4];
periCueOnsetInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, periCueOnsetWindowOffset);

periArrayOnsetWindowOffset = [-0.4 0.4];
periArrayOnsetInd = getTimeLogicalWithTolerance(EL.arrayOnsetHoldBalLfp.t, periArrayOnsetWindowOffset);

periTargetDimWindowOffset = [-0.4 0.4];
periTargetDimInd = getTimeLogicalWithTolerance(EL.targetDimBalLfp.t, periTargetDimWindowOffset);

params.tapers = [2 3];
params.fpass = [5 60];
params.pad = 2;
params.Fs = D.lfpFs;
params.trialave = 1;
yBounds = params.fpass;
spikeEventAlignWindow = [0.4 0.4]; % should match others
cBounds = [0 0.16];
cDiffBounds = [-0.1 0.1];

channelInd = nChannels+1; % CAR

inRFLoc = 3;
exRFLoc = 1;

cueResponses = cell(nLoc, 1);
arrayOnsetHoldResponses = cell(nLoc, 1);
targetDimResponses = cell(nLoc, 1);

for m = 1:nLoc
    cueOnsetLfpCurrent = squeeze(EL.cueOnsetLfp.lfp(channelInd,UE.cueLoc == m,:))'; % each column is a trial
    if isempty(cueOnsetLfpCurrent)
        continue;
    end

    arrayOnsetHoldLfpCurrent = squeeze(EL.arrayOnsetHoldBalLfp.lfp(channelInd,UE.cueLocHoldBal == m,:))'; % each column is a trial
    targetDimLfpCurrent = squeeze(EL.targetDimBalLfp.lfp(channelInd,UE.cueLocHoldBal == m,:))'; % each column is a trial

    cueResponses{m} = cueOnsetLfpCurrent(periCueOnsetInd,:);
    arrayOnsetHoldResponses{m} = arrayOnsetHoldLfpCurrent(periArrayOnsetInd,:);
    targetDimResponses{m} = targetDimLfpCurrent(periTargetDimInd,:);
end

nUnits = numel(D.allMUAStructs);
CDiffCueOnset = nan(nUnits, 13, 56);
CDiffArrayOnsetHold = nan(nUnits, 13, 56);
CDiffTargetDim = nan(nUnits, 13, 56);

for unitInd = 1:nUnits
    unitIDChar = D.allMUAStructs{unitInd}.unitIDChar;
    spikeTs = D.allMUAStructs{unitInd}.ts;
    
    figure_tr_inch(16, 10);
    
    suptitle(sprintf('SFC Unit %s - Channel %d\n', unitIDChar, channelInd));
        
    subaxis(3, 3, 1);
    alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.cueOnsetByLoc{inRFLoc}, spikeEventAlignWindow);
    [C,phi,S12,S1,S2,t,f] = cohgramcpt(cueResponses{inRFLoc}, alignedSpikeTs, movingWin, params);
    imagesc(t + periCueOnsetWindowOffset(1), f, C');
    set(gca, 'YDir', 'normal');
    xlim(periCueOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
    ylim(yBounds);
    caxis(cBounds);
    colorbar;
    CInRF = C;
    
    
    subaxis(3, 3, 4);
    alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.cueOnsetByLoc{exRFLoc}, spikeEventAlignWindow);
    [C,phi,S12,S1,S2,t,f] = cohgramcpt(cueResponses{exRFLoc}, alignedSpikeTs, movingWin, params);
    imagesc(t + periCueOnsetWindowOffset(1), f, C');
    set(gca, 'YDir', 'normal');
    xlim(periCueOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
    ylim(yBounds);
    caxis(cBounds);
    colorbar;
    CExRF = C;

    subaxis(3, 3, 7);
    imagesc(t + periCueOnsetWindowOffset(1), f, CInRF' - CExRF');
    set(gca, 'YDir', 'normal');
    xlim(periCueOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
    ylim(yBounds);
    caxis(cDiffBounds);
    colormap(gca, getCoolWarmMap());
    colorbar;
    CDiffCueOnset(unitInd,:,:) = CInRF - CExRF;
    
    subaxis(3, 3, 2);
    alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.arrayOnsetHoldBalByLoc{inRFLoc}, spikeEventAlignWindow);
    [C,phi,S12,S1,S2,t,f] = cohgramcpt(arrayOnsetHoldResponses{inRFLoc}, alignedSpikeTs, movingWin, params);
    imagesc(t + periArrayOnsetWindowOffset(1), f, C');
    set(gca, 'YDir', 'normal');
    xlim(periArrayOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
    ylim(yBounds);
    caxis(cBounds);
    colorbar;
    CInRF = C;
    
    subaxis(3, 3, 5);
    alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.arrayOnsetHoldBalByLoc{exRFLoc}, spikeEventAlignWindow);
    [C,phi,S12,S1,S2,t,f] = cohgramcpt(arrayOnsetHoldResponses{exRFLoc}, alignedSpikeTs, movingWin, params);
    imagesc(t + periArrayOnsetWindowOffset(1), f, C');
    set(gca, 'YDir', 'normal');
    xlim(periArrayOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
    ylim(yBounds);
    caxis(cBounds);
    colorbar;
    CExRF = C;

    subaxis(3, 3, 8);
    imagesc(t + periArrayOnsetWindowOffset(1), f, CInRF' - CExRF');
    set(gca, 'YDir', 'normal');
    xlim(periArrayOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
    ylim(yBounds);
    caxis(cDiffBounds);
    colormap(gca, getCoolWarmMap());
    colorbar;
    CDiffArrayOnsetHold(unitInd,:,:) = CInRF - CExRF;
    
    subaxis(3, 3, 3);
    alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.targetDimBalByLoc{inRFLoc}, spikeEventAlignWindow);
    [C,phi,S12,S1,S2,t,f] = cohgramcpt(targetDimResponses{inRFLoc}, alignedSpikeTs, movingWin, params);
    imagesc(t + periTargetDimWindowOffset(1), f, C');
    set(gca, 'YDir', 'normal');
    xlim(periTargetDimWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
    ylim(yBounds);
    caxis(cBounds);
    colorbar;
    CInRF = C;
    
    subaxis(3, 3, 6);
    alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.targetDimBalByLoc{exRFLoc}, spikeEventAlignWindow);
    [C,phi,S12,S1,S2,t,f] = cohgramcpt(targetDimResponses{exRFLoc}, alignedSpikeTs, movingWin, params);
    imagesc(t + periTargetDimWindowOffset(1), f, C');
    set(gca, 'YDir', 'normal');
    xlim(periTargetDimWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
    ylim(yBounds);
    caxis(cBounds);
    colorbar;
    CExRF = C;

    subaxis(3, 3, 9);
    imagesc(t + periTargetDimWindowOffset(1), f, CInRF' - CExRF');
    set(gca, 'YDir', 'normal');
    xlim(periTargetDimWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
    ylim(yBounds);
    caxis(cDiffBounds);
    colormap(gca, getCoolWarmMap());
    colorbar;
    CDiffTargetDim(unitInd,:,:) = CInRF - CExRF;
    
    drawnow;
    
    plotFileName = sprintf('%s/%s-FP%03d-%s-sfc-v%d.png', ...
            processedDataDir, fileNamePrefix, channelInd, unitIDChar, v);
    if channelInd == nChannels + 1
        plotFileName = sprintf('%s/%s-FPCAR-%s-sfc-v%d.png', ...
                processedDataDir, fileNamePrefix, unitIDChar, v);
    end
    fprintf('Saving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end

%% plot mean difference SFC InRF-ExRF

cDiffBounds = [-0.02 0.02];

figure_tr_inch(15, 5);

subaxis(1, 3, 1);
imagesc(t + periCueOnsetWindowOffset(1), f, squeeze(mean(CDiffCueOnset, 1))');
set(gca, 'YDir', 'normal');
xlim(periCueOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
ylim(yBounds);
caxis(cDiffBounds);
colormap(gca, getCoolWarmMap());

subaxis(1, 3, 2);
imagesc(t + periArrayOnsetWindowOffset(1), f, squeeze(mean(CDiffArrayOnsetHold, 1))');
set(gca, 'YDir', 'normal');
xlim(periArrayOnsetWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
ylim(yBounds);
caxis(cDiffBounds);
colormap(gca, getCoolWarmMap());

subaxis(1, 3, 3);
imagesc(t + periTargetDimWindowOffset(1), f, squeeze(mean(CDiffTargetDim, 1))');
set(gca, 'YDir', 'normal');
xlim(periTargetDimWindowOffset + [movingWin(1)/2 -movingWin(1)/2]);
ylim(yBounds);
caxis(cDiffBounds);
colormap(gca, getCoolWarmMap());

plotFileName = sprintf('%s/%s-FP%03d-allMUA-meanSfcDiff-v%d.png', ...
        processedDataDir, fileNamePrefix, channelInd, v);
if channelInd == nChannels + 1
    plotFileName = sprintf('%s/%s-FPCAR-allMUA-meanSfcDiff-v%d.png', ...
            processedDataDir, fileNamePrefix, v);
end
fprintf('Saving figure to file %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%%
return;
% below code has not been updated after the latest refactoring of EL vars

















%% image plot, bipolar reference

xBounds = [-0.4 0.4];
yBounds = [1 nChannels-1];
cBounds = [-0.02 0.02];

baselineWindowOffset = [-0.3 0];
baselineInd = getTimeLogicalWithTolerance(EL.cueOnsetT, baselineWindowOffset);

figure_tr_inch(20, 10); 
clf;
set(gcf, 'renderer', 'painters');

for m = 1:nLoc
    
preCueBaseline = diff(squeeze(mean(mean(EL.cueOnsetLfp(:,UE.cueLoc == m,baselineInd), 2), 3)));
    
subaxis(nLoc, 5, (m-1)*5 + 1, 'SH', 0.02, 'MT', 0.06);
hold on;
imagesc(EL.enterFixationT, 1:nChannels-1, diff(squeeze(mean(EL.enterFixationLfp(:,UE.cueLoc == m,:), 2))) - preCueBaseline);
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
set(gca, 'YDir', 'reverse');
xlim(xBounds);
ylim(yBounds);
caxis(cBounds);
set(gca, 'box', 'off');
if m == 1
    title('Start Fixation');
end
if m == nLoc
    xlabel('Time from Start Fixation (s)');
end
ylabel('Channel Number');

subaxis(nLoc, 5, (m-1)*5 + 2);
hold on;
imagesc(EL.cueOnsetT, 1:nChannels-1, diff(squeeze(mean(EL.cueOnsetLfp(:,UE.cueLoc == m,:), 2))) - preCueBaseline);
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
set(gca, 'YDir', 'reverse');
xlim(xBounds);
ylim(yBounds);
caxis(cBounds);
set(gca, 'box', 'off');
if m == 1
    title('Cue Onset');
end
if m == nLoc
    xlabel('Time from Cue Onset (s)');
end

subaxis(nLoc, 5, (m-1)*5 + 3);
hold on;
imagesc(EL.arrayOnsetT, 1:nChannels-1, diff(squeeze(mean(EL.arrayOnsetHoldLfp(:,UE.cueLocHold == m,:), 2))) - preCueBaseline);
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
set(gca, 'YDir', 'reverse');
xlim(xBounds);
ylim(yBounds);
caxis(cBounds);
set(gca, 'box', 'off');
if m == 1
    title('Array Onset');
end
if m == nLoc
    xlabel('Time from Array Onset (s)');
end

subaxis(nLoc, 5, (m-1)*5 + 4);
hold on;
imagesc(EL.targetDimT, 1:nChannels-1, diff(squeeze(mean(EL.targetDimLfp(:,UE.cueLocHold == m,:), 2))) - preCueBaseline);
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
set(gca, 'YDir', 'reverse');
xlim(xBounds);
caxis(cBounds);
ylim(yBounds);
set(gca, 'box', 'off');
if m == 1
    title('Target Dimming');
end
if m == nLoc
    xlabel('Time from Target Dimming (s)');
end

subaxis(nLoc, 5, (m-1)*5 + 5);
hold on;
imagesc(EL.exitFixationT, 1:nChannels-1, diff(squeeze(mean(EL.exitFixationLfp(:,UE.cueLoc == m,:), 2))) - preCueBaseline);
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
set(gca, 'YDir', 'reverse');
xlim(xBounds);
ylim(yBounds);
caxis(cBounds);
set(gca, 'box', 'off');
if m == 1
    title('Exit Fixation');
end
if m == nLoc
    xlabel('Time from Exit Fixation (s)');
end

end % end for each location

axBig = axes('Position', [0.04 0.035 0.92 0.93], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on');
title(axBig, sprintf('Session %s', sessionName), 'FontSize', 14);

plotFileName = sprintf('%s/%s-allFP-evokedLfps-%s-sepBipolarColorPlot-v%d.png', ...
        processedDataDir, blockName, ref, v);
export_fig(plotFileName, '-nocrop');

%% CSD
meanCueOnsetLfpByLoc = cell(nLoc, 1);
preCueBaselineIndices = getTimeLogicalWithTolerance(EL.cueOnsetT, EL.preCueBaselineWindowOffset);

yVals = 2:nChannels-1;
yBounds = yVals([1 end]) + [-0.5 0.5];
cBounds = [-0.02 0.02];
for m = 1:nLoc
    if all(isnan(meanCueOnsetLfpByLoc{m}(:)))
        continue;
    end
    meanCueOnsetLfpByLoc{m} = squeeze(mean(EL.cueOnsetLfp(:,UE.cueLoc == m,:), 2));
    cueOnsetCSD = nan(numel(yVals), size(meanCueOnsetLfpByLoc{m}, 2));
    figure;
    hold on;
    for j = 1:numel(yVals)
        ji = yVals(j);
        cueOnsetCSD(j,:) = meanCueOnsetLfpByLoc{m}(ji+1,:) - 2 * meanCueOnsetLfpByLoc{m}(ji,:) + meanCueOnsetLfpByLoc{m}(ji-1,:);
        cueOnsetCSD(j,:) = cueOnsetCSD(j,:) - mean(cueOnsetCSD(j,preCueBaselineIndices));
    end
    imagesc(EL.cueOnsetT, yVals, cueOnsetCSD);
    plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
    set(gca, 'YDir', 'reverse');
    xlim(xBounds);
    ylim(yBounds);
    caxis(cBounds);
    title(sprintf('P%d', m));
end

stop

%% restructure EL.cueOnsetLfpByLoc -> cueOnsetLfpByLoc2
cueOnsetLfpByLoc2 = cell(nLoc, 1);
nTime = numel(EL.cueOnsetT);
for j = 1:nChannels
    saveFileName = sprintf('%s/%s-%s-evokedLfp-v%d.mat', ...
            processedDataDir, D.lfpNames{j}, blockName, v);
    EL = load(saveFileName);
    for m = 1:nLoc
        nTrial = size(EL.cueOnsetLfpByLoc{m}, 1);
        if isempty(cueOnsetLfpByLoc2{m})
            cueOnsetLfpByLoc2{m} = nan(nChannels, nTrial, nTime);
        end
        for i = 1:nTrial
            cueOnsetLfpByLoc2{m}(j,:,:) = EL.cueOnsetLfpByLoc{m};
        end
    end
end

%%
specgramSlidingWindow = [0.4 0.05];
params.tapers = [1 1];
params.fpass = [4 30];
params.pad = 1;
params.trialave = 1;

j = 10;
saveFileName = sprintf('%s/%s-%s-evokedLfp-v%d.mat', ...
        processedDataDir, D.lfpNames{j}, blockName, v);
EL = load(saveFileName);

for m = 1:nLoc
    [S,t,f] = mtspecgramc(EL.arrayOnsetHoldLfpByLoc{m}', specgramSlidingWindow, params);
    figure;
    imagesc(t + EL.arrayOnsetWindowOffset(1), f, 10*log10(S'));
    set(gca, 'YDir', 'normal');
    xlim([-0.5 0.5]);
    colorbar;
    caxis([-62 -59]);
end

%% set up fieldtrip
cfg = struct();
cfg.padding = 0;
cfg.demean = 'yes';
cfg.baselinewindow = [-0.3 0];
cfg.detrend = 'yes';
cfg.reref = 'yes';
cfg.refchannel = 'all';
cfg.refmethod = 'avg';
cueOnsetLfp = convertToFieldTrip(D, adjLfpsClean, UE);
cueOnsetLfpPP = ft_preprocessing(cfg, cueOnsetLfp);

%%
cfg = struct();
cfg.method = 'tfr';
cfg.output = 'pow';
cfg.foilim = [5 80];
cfg.toi = -0.7:0.05:0.7;
tfrByLoc = cell(nLoc, 1);
tic;
for m = 1:nLoc
    fprintf('P%d:\n', m);
    cfg.trials = UE.cueLoc == m;
    tfrByLoc{m} = ft_freqanalysis(cfg, cueOnsetLfpPP);
end
toc

saveFileName = sprintf('%s/%s-allFP-evokedLfp-cueOnsetTFR-v%d.mat', ...
        processedDataDir, blockName, v);
save(saveFileName, 'tfrByLoc');

%%

% for i = 1:nChannels
%     figure;
%     imagesc(tfr.time, tfr.freq, squeeze(tfr.powspctrm(i,:,:)));
%     set(gca, 'YDir', 'normal');
%     colorbar;
% end

%%
load(saveFileName);

cfg              = [];
cfg.baseline     = [-0.3 0]; 
cfg.baselinetype = 'relchange'; 
cfg.maskstyle    = 'saturation';	
cfg.zlim         = [-1 3];	
cfg.interactive  = 'no';
cfg.ylim         = [5 30];
for i = 1:nChannels
    figure_tr_inch(10, 10);
    cfg.channel = i;
    for m = 1:nLoc
        subaxis(2, 2, m);
        ft_singleplotTFR(cfg, tfrByLoc{m});
        title(sprintf('P%d', m));
    end
    suptitle(sprintf('Channel %d', i));
end


%%
cfg = struct();
cfg.toi = -0.7:0.05:0.7;
cfg.vartrllength = 1;
erByLoc = cell(nLoc, 1);
tic;
for m = 1:nLoc
    fprintf('P%d:\n', m);
    cfg.trials = UE.cueLoc == m;
    erByLoc{m} = ft_timelockanalysis(cfg, cueOnsetLfpPP);
end
toc

saveFileName = sprintf('%s/%s-allFP-evokedLfp-cueOnsetER-v%d.mat', ...
        processedDataDir, blockName, v);
save(saveFileName, 'erByLoc');

%% fieldtrip: plot evoked response
load(saveFileName);

cfg              = [];
cfg.baseline     = [-0.3 0]; 
cfg.baselinetype = 'absolute'; 
cfg.maskstyle    = 'saturation';	
cfg.interactive  = 'no';
cfg.xlim         = [-0.3 0.3]; 
cfg.ylim         = [-0.02 0.02];
for i = 1:nChannels
    figure_tr_inch(10, 10);
    cfg.channel = i;
    for m = 1:nLoc
        subaxis(2, 2, m);
        cfg.graphcolor = cols(m,:);
        ft_singleplotER(cfg, erByLoc{m});
        hold on;
        plot([0 0], cfg.ylim, 'Color', 0.3*ones(3, 1));
        title(sprintf('P%d', m), 'FontSize', 14);
        
    end
    suptitle(sprintf('Channel %d', i));
end

% results not consistent with above

%%


%% try motif discovery

channelInd = 25;
m = 3;
postCueOnsetInd = 901:1250;

dataOrig = squeeze(EL.cueOnsetLfp(channelInd,UE.cueLoc == m,:))'; % each column is a trial
hiCutoffFreqCustom = 80; % low-pass filter at 80 Hz
bFirLowPassCustom = fir1(3*fix(Fs/hiCutoffFreqCustom), hiCutoffFreqCustom/(Fs/2), 'low');
dataFilt = filtfilt(bFirLowPassCustom, 1, dataOrig);
data = dataFilt(postCueOnsetInd,:); % concatenate trials, column append to column
data = data(:);

%%
tic; 
segmentSize = 100;
[matrixProfile, profileIndex, motifIdxs] = interactiveMatrixProfileVer2(data, segmentSize); 
toc

% relative time
mod(sort(motifIdxs{1,2}) + segmentSize - 1, numel(postCueOnsetInd))
mod(sort(motifIdxs{2,2}) + segmentSize - 1, numel(postCueOnsetInd))
mod(sort(motifIdxs{3,2}) + segmentSize - 1, numel(postCueOnsetInd))

%%
yBounds = [-1 1] * 0.05;
figure;
hold on;
plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
j = channelInd;
plot(EL.cueOnsetT, dataFilt(:,1:3), 'Color', cols(m,:));
xlim(EL.cueOnsetT([1 end]));
ylim(yBounds);
set(gca, 'box', 'off');
xlabel('Time from Cue Onset (s)');

%% plot after filtering
postCueResponseOffset = [0.2 0.5];
postCueResponseInd = getTimeLogicalWithTolerance(EL.cueOnsetT, postCueResponseOffset);

baselineWindowOffset = [-0.3 0];
baselineInd = getTimeLogicalWithTolerance(EL.cueOnsetT, baselineWindowOffset);

for i = 1:nChannels
    channelInd = i;
    cueOnsetLfpFilt = cell(nLoc, 1);
    preCueBaselineMeanResponse = nan(nLoc, 1);
    for m = 1:nLoc
        dataOrig = squeeze(EL.cueOnsetLfp(channelInd,UE.cueLoc == m,:))'; % each column is a trial
        hiCutoffFreqCustom = 80; % low-pass filter at 80 Hz
        bFirLowPassCustom = fir1(3*fix(Fs/hiCutoffFreqCustom), hiCutoffFreqCustom/(Fs/2), 'low');
        cueOnsetLfpFilt{m} = filtfilt(bFirLowPassCustom, 1, dataOrig);
        preCueBaselineMeanResponse(m) = mean(mean(cueOnsetLfpFilt{m}(baselineInd,:), 1), 2);
    end

    yBounds = [-1 1] * 0.015;
    figure_tr_inch(12, 6);
    hold on;
    plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
    j = channelInd;
    for m = 1:nLoc
        plot(EL.cueOnsetT, mean(cueOnsetLfpFilt{m}, 2) - preCueBaselineMeanResponse(m), 'Color', cols(m,:), 'LineWidth', 2);
    end
    xlim(EL.cueOnsetT([1 end]));
    ylim(yBounds);
    set(gca, 'box', 'off');
    xlabel('Time from Cue Onset (s)');
    title(sprintf('Channel %d', i));
end

%%
params.tapers = [2 3];
params.fpass = [5 80];
params.pad = 2;
specgramMovingWindow = [0.2 0.05];
[S,t,f,Serr] = mtspecgramc(cueOnsetLfpFilt{3}, specgramMovingWindow, params);

%%
figure;
imagesc(t - 0.7, f, 10*log10(S'));
set(gca, 'YDir', 'normal');

