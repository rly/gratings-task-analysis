function lfpVEPMapping(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, channelsToLoad, ref)
% LFP VEP Mapping, all channels on a probe
% can't really do one channel at a time because of Common Average
% Referencing

%% LFP VEP mapping parameters
% ref = 'RAW'; % CAR, BIP, or some other code, e.g. RAW. CAR recommended

%% setup and load data
v = 11;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('VEP Mapping Analysis - LFP\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('LFP Channel to Load: %d\n', channelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Reference: %s\n', ref);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

%% input check
% assert(numel(channelsToLoad) == 1);
assert(numel(channelsToLoad) > 1);
assert(numel(channelsToLoad) <= 32);

%% load recording information
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, channelsToLoad, 'VEPM', 'LFP_VEPM', 1, 1);
sessionName = R.sessionName;
areaName = R.areaName;

%% extract flash events
% Flash mapping events pre 20170202
% EVT05 - begin pre-cue period
% EVT06 - flash onset
% EVT07 - FP dim
% EVT08 - juice onset

% Flash mapping events 20170202 and on
% EVT02 - begin pre-cue period
% EVT03 - flash onset
% EVT07 - FP dim
% EVT08 - juice onset

if strcmp(sessionName, 'M20170127') || strcmp(sessionName, 'M20170130') || strcmp(sessionName, 'M20170201')
    origFlashEvents = D.events{6};
    preFlashesEvents = D.events{5};
else
    origFlashEvents = D.events{3};
    preFlashesEvents = D.events{2};
end

%%
periFlashWindowOffset = [-0.25 0.3]; % seconds around flash
baselineWindowOffset = [-0.2 0]; % seconds - 200ms minimum between flashes and before first flash
expandedPlotWindowOffset = [-0.5 1.5]; % seconds - 150ms minimum between flashes and before first flash
% NOTE: for some early sessions, there were only 3 flashes per trial and so
% expandedPlot should only go until 1.5 s

plotFileNamePrefix = sprintf('%s_%s_ch%d-ch%d_%s-FPall', sessionName, areaName, channelsToLoad([1 end]), blockName);

%% preprocess LFPs
Fs = D.lfpFs;
nChannels = D.nLfpCh;

D.adjLfpsClean = interpolateLfpOverSpikeTimes(D.adjLfps, channelsToLoad, Fs, D.allMUAStructs);

% TODO use isEventOutlier
hiCutoffFreq = 100;
[channelDataNorm,flashEventsClean,isEventOutlier,isNoisyChannel] = preprocessLfps(D.adjLfpsClean, Fs, D.lfpNames, origFlashEvents, ...
        processedDataDir, blockName, hiCutoffFreq, 1, v);
% channelDataNorm = D.adjLfpsClean;
% flashEventsClean = origFlashEvents;
D.adjLfps = [];
D.adjLfpsClean = [];
nFlashes = numel(flashEventsClean);

%% calculate and plot per-trial baseline
fPlot = figure_tr_inch(6, 10);
[tBaseline,averageBaselineResponses] = computeResponsesInWindow(channelDataNorm, preFlashesEvents, ...
        baselineWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaseline, averageBaselineResponses, 'yScale', 500);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');

plotFileName = sprintf('%s/%s_baseline.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName, '-nocrop');

%% plot baseline responses overlaid
figure_tr_inch(6, 5, 6, 0);
hold on;
plot(tBaseline, averageBaselineResponses - repmat(mean(averageBaselineResponses, 2), 1, size(averageBaselineResponses, 2)));
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');

plotFileName = sprintf('%s/%s_baselineOverlaid.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName, '-nocrop');

%% plot expanded range around trial onset
fPlot = figure_tr_inch(16, 5, 0, 6);
[tBaselineExp,averageBaselineResponsesExp] = computeResponsesInWindow(channelDataNorm, preFlashesEvents, ...
        expandedPlotWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaselineExp, averageBaselineResponsesExp, 'yScale', 100);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');
xlim(expandedPlotWindowOffset);

plotFileName = sprintf('%s/%s_baselineExpanded.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName, '-nocrop');

%% plot flash responses overlaid
figure_tr_inch(6, 5, 6, 0);
hold on;
plot(tBaselineExp, averageBaselineResponsesExp - repmat(mean(averageBaselineResponsesExp, 2), 1, size(averageBaselineResponsesExp, 2)));
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');

plotFileName = sprintf('%s/%s_baselineExpandedOverlaid.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName, '-nocrop');

%% plot expanded range around trial onset after cross-channel CAR
fPlot = figure_tr_inch(16, 5, 0, 6);
[tBaselineExp,averageBaselineResponsesCARExp] = computeResponsesInWindow(channelDataNorm - nanmean(channelDataNorm, 1), preFlashesEvents, ...
        expandedPlotWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaselineExp, averageBaselineResponsesCARExp, 'yScale', 500);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel (CAR)');
xlim(expandedPlotWindowOffset);

plotFileName = sprintf('%s/%s_baselineCARExpanded.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName, '-nocrop');

%% plot expanded range around trial onset after cross-channel BIP
% works only for one channel right now
fPlot = figure_tr_inch(16, 5, 0, 6);
[tBaselineExp,averageBaselineResponsesBIPExp] = computeResponsesInWindow(channelDataNorm(2:end,:) - channelDataNorm(1:end-1,:), preFlashesEvents, ...
        expandedPlotWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaselineExp, averageBaselineResponsesBIPExp, 'yScale', 100);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel (BIP)');
xlim(expandedPlotWindowOffset);

plotFileName = sprintf('%s/%s_baselineBIPExpanded.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName, '-nocrop');

%% plot per channel responses after outlier removal
fprintf('\nProcessing channel responses after outlier removal...\n');

% recompute baseline (?)

startIndices = round((flashEventsClean + periFlashWindowOffset(1)) * Fs); 
endIndices = startIndices + round(diff(periFlashWindowOffset) * Fs) - 1;
t = periFlashWindowOffset(1):1/Fs:periFlashWindowOffset(2)-1/Fs;
nTime = numel(t);
assert(all(nTime == (endIndices - startIndices + 1)));

if strcmp(ref, 'BIP') % not the best method -- may be better to run this script on each probe separately
    if nChannels == 32
        nChannels = 31;
    end
end
meanLfpAcrossChannels = nanmean(channelDataNorm, 1);
responses = nan(nChannels, nTime, nFlashes);
for j = 1:nChannels
    % normalize each channel by mean baseline across trials and time and
    % std across time of mean baseline across trials
    % TODO: bar plot of average baseline response and std by channel
    fprintf('Processing %s...\n', D.lfpNames{j});
    
    % subtracting baseline not yet fully implemented
%     channelDataBaselined = (channelData - nanmean(averageBaselineResponses(j,:))) / nanstd(averageBaselineResponses(j,:));

    if strcmp(ref, 'CAR')
        channelData = channelDataNorm(j,:) - meanLfpAcrossChannels;
    elseif strcmp(ref, 'BIP')
        channelData = channelDataNorm(j+1,:) - channelDataNorm(j,:);
    else % RAW
        channelData = channelDataNorm(j,:);
    end
    channelDataRenorm = (channelData - nanmean(channelData)) / nanstd(channelData);
    
    % can vectorize???
    for i = 1:nFlashes
        responses(j,:,i) = channelDataRenorm(startIndices(i):endIndices(i));
    end
end

clear channelData channelDataNorm;

%% subtract out baseline
% units are standard deviations from baselined mean
averageResponse = mean(responses, 3); % average across flashes

% subtract mean baseline activity (-0.1, 0] seconds before flash
flashBaselineWindowOffset = [-0.1 0];

% TODO make so that baseline can have arbitrary end time
indexFlashTime = -round(periFlashWindowOffset(1) * Fs);
indexStartBaselineFlashTime = -round((periFlashWindowOffset(1) - flashBaselineWindowOffset(1)) * Fs) + 1;
for j = 1:nChannels
    averageResponse(j,:) = averageResponse(j,:) - mean(averageResponse(j,indexStartBaselineFlashTime:indexFlashTime));
end

%% staggered channel line plot of average visually evoked LFP
xBounds = [-0.1 0.3];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(averageResponse)));
figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset) - 1]);
    end
end
xlim(periFlashWindowOffset);
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 12);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 12);
xlabel('Time from Flash Onset (s)');
title(sprintf('%s %s - LFP Response to Full-Field Flash (N = %d) (Ref: %s)', sessionName, areaName, nFlashes, ref));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, minYLim+1, '35-45 ms', 'Color', 'm');

plotFileName = sprintf('%s/%s_%s_lfpLines.png', processedDataDir, plotFileNamePrefix, ref);
export_fig(plotFileName, '-nocrop');

%% staggered channel color plot of average visually evoked LFP
figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.08);
hold on;
imagesc(t, 1:nChannels, averageResponse);
set(gca, 'YDir', 'reverse');
plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim([0.5 nChannels+0.5]);
xlabel('Time from Flash Onset (s)');
ylabel('Channel Number (1 = topmost)');
title(sprintf('%s %s - LFP Response to Full-Field Flash (N = %d) (Ref: %s)', sessionName, areaName, nFlashes, ref));
maxCAxis = max(abs(caxis));%*3/4;
caxis([-maxCAxis maxCAxis]);
colormap(getCoolWarmMap());
colorbar;
set(gca, 'FontSize', 12);

if (nChannels > 31 && strcmp(ref, 'BIP')) || (nChannels > 32 && strcmp(ref, 'CAR'))
    plot(xlim(), [nChannels/2+0.5 nChannels/2+0.5], 'Color', 0.3*ones(3, 1));
end

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, nChannels, '35-45 ms', 'Color', 'm');

plotFileName = sprintf('%s/%s_%s_lfpColor.png', processedDataDir, plotFileNamePrefix, ref);
export_fig(plotFileName, '-nocrop');

%% plot boxes around areas abs() > thresh over last plot and re-save
boxAbsThresh = 0.25;
strongResponseGroups = bwlabel(abs(averageResponse) > boxAbsThresh, 4);
strongResponseGroupBoundaries = bwboundaries(strongResponseGroups);
for k = 1:length(strongResponseGroupBoundaries)
   boundary = strongResponseGroupBoundaries{k};
   plot(t(boundary(:,2)), boundary(:,1), 'Color', 0.3*ones(3, 1), 'LineWidth', 2)
   
   minBY = min(boundary(:,1));
   maxBY = max(boundary(:,1));
   text(-0.01, minBY, sprintf('%d', minBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
   text(-0.01, maxBY, sprintf('%d', maxBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
end

plotFileName = sprintf('%s/%s_%s_lfpColorBounds.png', processedDataDir, plotFileNamePrefix, ref);
export_fig(plotFileName, '-nocrop');

%% plot vertical line plot of mean activity between 35 and 45 ms after flash
% TODO mark which ones are more than 2 SDs from baseline for that channel

earlyActivityWindowOffset = [0.035 0.045];
earlyActivityTLogical = getTimeLogicalWithTolerance(t, earlyActivityWindowOffset);
meanEarlyActivity = mean(averageResponse(:,earlyActivityTLogical), 2);
channelIndices = 1:nChannels;

figure_tr_inch(7, 8);
subaxis(1, 1, 1, 'ML', 0.12, 'MB', 0.11, 'MR', 0.05);
hold on;
plot(meanEarlyActivity, channelIndices, '.--', 'MarkerSize', 25, 'LineWidth', 2, 'Color', lines(1));
plot([0 0], channelIndices([1 end]) + [-1 1], '-', 'Color', 0.3*ones(3, 1));
set(gca, 'FontSize', 16);
set(gca, 'YDir', 'reverse');
xlabel('Mean Early Response to Flash (Overall SDs)');
ylabel('Channel Number');
grid on;
ylim(channelIndices([1 end]) + [-1 1]);
xlim([-1 1] * max(abs(xlim())));
title(sprintf('%s %s - Mean Early Response (%d-%d ms)', sessionName, areaName, earlyActivityWindowOffset * 1000));

plotFileName = sprintf('%s/%s_%s_meanEarlyResponse.png', processedDataDir, plotFileNamePrefix, ref);
export_fig(plotFileName, '-nocrop');
