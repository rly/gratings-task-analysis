% function lfpVEPMapping(pl2FileName, sessionName, areaName, blockNames, blockInds)

clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';

sessionInd = 8;

%% load recording information
recordingInfo = readRecordingInfo();
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% LFP VEP mapping parameters
ref = 'RAW'; % CAR, BIP, or some other code, e.g. RAW. CAR recommended
doOutlierCheckPlot = 0;

%% setup and load data
fprintf('\n-------------------------------------------------------\n');
fprintf('VEP Mapping Analysis - LFPs\n');
fprintf('Loading %s...\n', pl2FilePath);

tic;
isLoadSpikes = 1;
isLoadLfp = 1;
isLoadSpkc = 0;
isLoadDirect = 0;
D = loadPL2(pl2FilePath, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad); 

processedDataDir = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

blockName = strjoin(blockNames(vepmIndices), '-');

fprintf('Reference: %s\n', ref);

%% remove spike and event times not during RFM task to save memory
D = trimSpikeTimesAndEvents(D, vepmIndices);
% TODO should trim LFPs too!!!
% TODO show individual trials

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
nFlashes = numel(origFlashEvents);

%%
periFlashWindowOffset = [-0.25 0.3]; % ms around flash
baselineWindowOffset = [-0.2 0]; % ms - 200ms minimum between flashes and before first flash
expandedPlotWindowOffset = [-0.5 1.8]; % ms - 150ms minimum between flashes and before first flash

maxAbsYNonOutlier = 1000;
outlierMaxSD = 6;
outlierCheckWindowOffset = periFlashWindowOffset; 

Fs = D.lfpFs;
nChannels = D.nLfpCh;

plotFileNamePrefix = sprintf('%s_%s_ch%d-ch%d_%s-FPall', sessionName, areaName, lfpChannelsToLoad([1 end]), blockName);

%% low-pass filter data even further - this is slow
% ideally don't low-pass at 300 and then again at another freq
tic;
fprintf('Low-pass filtering data...');
hiCutoffFreqCustom = 100; % low-pass filter at 100 Hz
bFirLowPassCustom = fir1(3*fix(Fs/hiCutoffFreqCustom), hiCutoffFreqCustom/(Fs/2), 'low');
D.adjLfps(isnan(D.adjLfps)) = 0;
D.adjLfpsLP = filtfilt(bFirLowPassCustom, 1, D.adjLfps');
D.adjLfpsLP = D.adjLfpsLP';
D.adjLfps = [];
fprintf('... done (%0.2f s).\n', toc);

%% calculate and plot per-trial baseline
fPlot = figure_tr_inch(6, 10);
[tBaseline,averageBaselineResponses] = computeResponsesInWindow(D.adjLfpsLP, preFlashesEvents, ...
        baselineWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaseline, averageBaselineResponses, 'yScale', 500);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');

plotFileName = sprintf('%s/%s_baseline.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName);

%% plot baseline responses overlaid
figure_tr_inch(6, 5, 6, 0);
hold on;
plot(tBaseline, averageBaselineResponses - repmat(mean(averageBaselineResponses, 2), 1, size(averageBaselineResponses, 2)));
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');

plotFileName = sprintf('%s/%s_baselineOverlaid.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName);

%% plot expanded range around trial onset
fPlot = figure_tr_inch(16, 5, 0, 6);
[tBaselineExp,averageBaselineResponsesExp] = computeResponsesInWindow(D.adjLfpsLP, preFlashesEvents, ...
        expandedPlotWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaselineExp, averageBaselineResponsesExp, 'yScale', 100);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');
xlim(expandedPlotWindowOffset);

plotFileName = sprintf('%s/%s_baselineExpanded.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName);

%% plot flash responses overlaid
figure_tr_inch(6, 5, 6, 0);
hold on;
plot(tBaselineExp, averageBaselineResponsesExp - repmat(mean(averageBaselineResponsesExp, 2), 1, size(averageBaselineResponsesExp, 2)));
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel');

plotFileName = sprintf('%s/%s_baselineExpandedOverlaid.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName);

%% plot expanded range around trial onset after cross-channel CAR
fPlot = figure_tr_inch(16, 5, 0, 6);
[tBaselineExp,averageBaselineResponsesCARExp] = computeResponsesInWindow(D.adjLfpsLP - nanmean(D.adjLfpsLP, 1), preFlashesEvents, ...
        expandedPlotWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaselineExp, averageBaselineResponsesCARExp, 'yScale', 500);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel (CAR)');
xlim(expandedPlotWindowOffset);

plotFileName = sprintf('%s/%s_baselineCARExpanded.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName);

%% plot expanded range around trial onset after cross-channel BIP
% works only for one channel right now
fPlot = figure_tr_inch(16, 5, 0, 6);
[tBaselineExp,averageBaselineResponsesBIPExp] = computeResponsesInWindow(D.adjLfpsLP(2:end,:) - D.adjLfpsLP(1:end-1,:), preFlashesEvents, ...
        expandedPlotWindowOffset, Fs);
plotPeriEventTracesLinearArray(fPlot, tBaselineExp, averageBaselineResponsesBIPExp, 'yScale', 100);
xlabel('Time from Pre-Flashes Marker (s)');
title('Mean Baseline Activity on Each Channel (BIP)');
xlim(expandedPlotWindowOffset);

plotFileName = sprintf('%s/%s_baselineBIPExpanded.png', processedDataDir, plotFileNamePrefix);
export_fig(plotFileName);

%% detect events with outlier activity (e.g. amp saturated)
fprintf('\nProcessing channel responses and removing outliers...\n');

outlierCheckStartIndices = round((origFlashEvents + outlierCheckWindowOffset(1)) * Fs); 
outlierCheckEndIndices = outlierCheckStartIndices + round(diff(outlierCheckWindowOffset) * Fs) - 1;
outlierCheckT = outlierCheckWindowOffset(1):1/Fs:outlierCheckWindowOffset(2)-1/Fs;
nOutlierCheckTime = numel(outlierCheckT);
assert(all(nOutlierCheckTime == (outlierCheckEndIndices - outlierCheckStartIndices + 1)));

isFlashEventOutlier = false(nFlashes, 1);
fprintf('Checking for flashes with outlier activity (e.g. amp saturated)...');

% do CAR for outlier removal
meanLfpAcrossChannels = nanmean(D.adjLfpsLP, 1);
for j = 1:nChannels
    fprintf('Processing %s...\n', D.lfpNames{j});
    
    channelData = D.adjLfpsLP(j,:) - meanLfpAcrossChannels;
    
    % normalize each channel by its mean and standard deviation
    % units are now standard deviations away from the mean
    channelDataNorm = (channelData - nanmean(channelData)) / nanstd(channelData);
    
    % plot per channel responses before outlier removal
    if doOutlierCheckPlot
        fHandles(j) = figure_tr_inch(12,6);
        suptitle(sprintf('Outlier Check on Channel %d', j));
        subaxis(1,2,1);
        hold on;
    end
    
    for i = 1:nFlashes
        dataInOutlierCheckPeriod = channelDataNorm(outlierCheckStartIndices(i):outlierCheckEndIndices(i));
        % mark outlier flashes
        if any(isnan(dataInOutlierCheckPeriod))
            isFlashEventOutlier(i) = 1;
            fprintf('\tDetected outlier flash (index %d) because of NaN (%d time points) in trial\n', ...
                    i, sum(isnan(dataInOutlierCheckPeriod)));
        elseif any(abs(dataInOutlierCheckPeriod) > maxAbsYNonOutlier)
            isFlashEventOutlier(i) = 1;
            fprintf('\tDetected outlier flash (index %d) because absolute value of voltage in trial is too high (%d time points > %0.1f)\n', ...
                    i, sum(abs(dataInOutlierCheckPeriod) > maxAbsYNonOutlier), maxAbsYNonOutlier);
        elseif any(abs(dataInOutlierCheckPeriod) > outlierMaxSD)
            isFlashEventOutlier(i) = 1;
            fprintf('\tDetected outlier flash (index %d) because voltage in trial is too extreme (%d time points > %d SDs)\n', ...
                    i, sum(abs(dataInOutlierCheckPeriod) > outlierMaxSD), outlierMaxSD);
        end
    
        if doOutlierCheckPlot
            plot(outlierCheckT, dataInOutlierCheckPeriod);
            xlim(outlierCheckWindowOffset);
        end
    end
    
    if doOutlierCheckPlot
        title('Original Responses');
        xlabel('Time from Flash Onset');
        ylabel('SDs from Overall Mean');
        origYLim = ylim();
        plot([0 0], [-1000 1000], '-', 'Color', [0.5 0.5 0.5]);
        ylim(origYLim);
    end
end
clear channelData channelDataNorm;

% remove outlier flashes
fprintf('\tDetected %d outlier flashes. Removing those flashes.\n', sum(isFlashEventOutlier));
flashEvents = origFlashEvents;
flashEvents(isFlashEventOutlier) = [];
nFlashes = numel(flashEvents);

%% plot per channel responses after outlier removal
fprintf('\nProcessing channel responses after outlier removal...\n');

% recompute baseline (?)

startIndices = round((flashEvents + periFlashWindowOffset(1)) * Fs); 
endIndices = startIndices + round(diff(periFlashWindowOffset) * Fs) - 1;
t = periFlashWindowOffset(1):1/Fs:periFlashWindowOffset(2)-1/Fs;
nTime = numel(t);
assert(all(nTime == (endIndices - startIndices + 1)));

if strcmp(ref, 'BIP') % not the best method -- may be better to run this script on each probe separately
    if nChannels == 32
        nChannels = 31;
    elseif nChannels == 64
        nChannels = 62;
    end
end
meanLfpAcrossChannels = nanmean(D.adjLfpsLP, 1);
responses = nan(nChannels, nTime, nFlashes);
for j = 1:nChannels
    % normalize each channel by mean baseline across trials and time and
    % std across time of mean baseline across trials
    % TODO: bar plot of average baseline response and std by channel
    fprintf('Processing %s...\n', D.lfpNames{j});
    
    % subtracting baseline not yet fully implemented
%     channelDataBaselined = (channelData - nanmean(averageBaselineResponses(j,:))) / nanstd(averageBaselineResponses(j,:));

    if strcmp(ref, 'CAR')
        channelData = D.adjLfpsLP(j,:) - meanLfpAcrossChannels;
    elseif strcmp(ref, 'BIP')
        if nChannels == 62 && j > 31 % separate the probes properly
            channelData = D.adjLfpsLP(j+2,:) - D.adjLfpsLP(j+1,:);
        else
            channelData = D.adjLfpsLP(j+1,:) - D.adjLfpsLP(j,:);
        end
    else
        channelData = D.adjLfpsLP(j,:);
    end
    channelDataNorm = (channelData - nanmean(channelData)) / nanstd(channelData);
    
    % can vectorize???
    for i = 1:nFlashes
        responses(j,:,i) = channelDataNorm(startIndices(i):endIndices(i));
    end
end

clear channelData channelDataNorm;

% plot for sanity check
if doOutlierCheckPlot
    for j = 1:nChannels
        figure(fHandles(j)); % recall the old figure handle
        subaxis(1,2,2);
        hold on;
        plot(t, squeeze(responses(j,:,:)));
        xlim(outlierCheckWindowOffset);
        title('After Outlier Removal');
        xlabel('Time from Flash Onset');
        origYLim = ylim();
        plot([0 0], [-1000 1000], '-', 'Color', [0.5 0.5 0.5]);
        ylim(origYLim);
    end
end

%% subtract out baseline
% units are standard deviations from baselined mean
averageResponse = mean(responses, 3);

flashBaselineWindowOffset = [-0.1 0];

% baseline correct based on (-0.1, 0] mean activity before flash
indexFlashTime = -round(periFlashWindowOffset(1) * Fs);
indexStartBaselineFlashTime = -round((periFlashWindowOffset(1) - flashBaselineWindowOffset(1)) * Fs) + 1;
for j = 1:nChannels
    averageResponse(j,:) = averageResponse(j,:) - mean(averageResponse(j,indexStartBaselineFlashTime:indexFlashTime));
end

%% staggered channel line plot of average visually evoked LFP
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
title(sprintf('Average LFP Response to Full-Field Flash (%d Flashes) (Ref: %s)', nFlashes, ref));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
ylim(origYLim);

% plot early latency line at 35ms
plot([0.035 0.035], origYLim, 'm-');
text(0.04, -(nChannels + 1), '35 ms', 'Color', 'm');

plotFileName = sprintf('%s/%s_%s_lfpLines.png', processedDataDir, plotFileNamePrefix, ref);
export_fig(plotFileName);

%% staggered channel color plot of average visually evoked LFP
figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.08);
hold on;
imagesc(t, 1:nChannels, averageResponse);
set(gca, 'YDir', 'reverse');
plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
xlim(periFlashWindowOffset);
ylim([0.5 nChannels+0.5]);
xlabel('Time from Flash Onset (s)');
ylabel('Channel Number (1 = topmost)');
title(sprintf('Average LFP Response to Full-Field Flash (%d Flashes) (Ref: %s)', nFlashes, ref));
maxCAxis = max(abs(caxis));%*3/4;
caxis([-maxCAxis maxCAxis]);
colormap(getCoolWarmMap());
colorbar;
set(gca, 'FontSize', 12);

if (nChannels > 31 && strcmp(ref, 'BIP')) || (nChannels > 32 && strcmp(ref, 'CAR'))
    plot(xlim(), [nChannels/2+0.5 nChannels/2+0.5], 'Color', 0.3*ones(3, 1));
end

% plot early latency line at 35ms
plot([0.035 0.035], [0 nChannels+1], 'm-');
text(0.04, nChannels+0.15, '35 ms', 'Color', 'm');

plotFileName = sprintf('%s/%s_%s_lfpColor.png', processedDataDir, plotFileNamePrefix, ref);
export_fig(plotFileName);

%% plot boxes around areas abs() > thresh over last plot and re-save
boxAbsThresh = 0.25;
if nChannels == 32 || nChannels == 31
    strongResponseGroups = bwlabel(abs(averageResponse) > boxAbsThresh, 4);
    strongResponseGroupBoundaries = bwboundaries(strongResponseGroups);
elseif nChannels == 64 % don't merge groups across probe boundaries
    strongResponseGroups1 = bwlabel(abs(averageResponse(1:32,:)) > boxAbsThresh, 4);
    strongResponseGroups2 = bwlabel(abs(averageResponse(33:64,:)) > boxAbsThresh, 4);
    strongResponseGroups2(strongResponseGroups2 > 0) = strongResponseGroups2(strongResponseGroups2 > 0) + max(strongResponseGroups1(:));
    strongResponseGroupBoundaries1 = bwboundaries([strongResponseGroups1; zeros(32, size(strongResponseGroups1, 2))]);
    strongResponseGroupBoundaries2 = bwboundaries([zeros(32, size(strongResponseGroups2, 2)); strongResponseGroups2]);
    strongResponseGroupBoundaries = [strongResponseGroupBoundaries1; strongResponseGroupBoundaries2];
elseif nChannels == 62 % don't merge groups across probe boundaries
    strongResponseGroups1 = bwlabel(abs(averageResponse(1:31,:)) > boxAbsThresh, 4);
    strongResponseGroups2 = bwlabel(abs(averageResponse(32:62,:)) > boxAbsThresh, 4);
    strongResponseGroups2(strongResponseGroups2 > 0) = strongResponseGroups2(strongResponseGroups2 > 0) + max(strongResponseGroups1(:));
    strongResponseGroupBoundaries1 = bwboundaries([strongResponseGroups1; zeros(31, size(strongResponseGroups1, 2))]);
    strongResponseGroupBoundaries2 = bwboundaries([zeros(31, size(strongResponseGroups2, 2)); strongResponseGroups2]);
    strongResponseGroupBoundaries = [strongResponseGroupBoundaries1; strongResponseGroupBoundaries2];
end
for k = 1:length(strongResponseGroupBoundaries)
   boundary = strongResponseGroupBoundaries{k};
   plot(t(boundary(:,2)), boundary(:,1), 'Color', 0.3*ones(3, 1), 'LineWidth', 2)
   
   minBY = min(boundary(:,1));
   maxBY = max(boundary(:,1));
   text(-0.01, minBY, sprintf('%d', minBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
   text(-0.01, maxBY, sprintf('%d', maxBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
end

plotFileName = sprintf('%s/%s_%s_lfpColorBounds.png', processedDataDir, plotFileNamePrefix, ref);
export_fig(plotFileName);



