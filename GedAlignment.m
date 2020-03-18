% CAR LFP data for VEP data

%'/Volumes/scratch/rly/gratings-task-analysis/processed_data//M20170324/LFP_VEPM/M20170324-ind11-PUL-ch1-ch32-vepm1-vepm2-vepm3-vepm4-vepm5-CAR-responses-v12.mat'
cd('/Volumes/scratch/rly/gratings-task-analysis/processed_data//M20170324/LFP_VEPM/')
load('M20170324-ind12-PUL-ch33-ch64-vepm1-vepm2-vepm3-vepm4-vepm5-CAR-responses-v12.mat')

processedDataRootDir = '/Volumes/scratch/rly/gratings-task-analysis/processed_data/';
dataDirRoot = '/Volumes/kastner/ryanly/McCartney/merged';
muaDataDirRoot = '/Volumes/scratch/rly/simple-mua-detection/processed_data/';
recordingInfoFileName = '/Users/labmanager/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
sessionInd = 12;
channelsToLoad = 33:64;
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, channelsToLoad, 'VEPM', 'LFP_VEPM', 0, 0, 1, 0);
sessionName = R.sessionName;
areaName = R.areaName;

if strcmp(sessionName, 'M20170127') || strcmp(sessionName, 'M20170130') || strcmp(sessionName, 'M20170201')
    origFlashEvents = D.events{6};
    preFlashesEvents = D.events{5};
else
    origFlashEvents = D.events{3};
    preFlashesEvents = D.events{2};
end

%[fixationAndLeverTimes,isMissingData] = getFixationAndLeverTimes(D, cueOnset, firstJuiceEvent, cueLoc, isHoldTrial, nLoc)

%%
% flash related activity 200ms windows
postFlashWindowOffset = [0.025 0.225]; % seconds after flash
baselineWindowOffset = [-.175 .025]; % seconds around preflashevent
saccadeWindowOffset = [-.425 -.225]; % fixation onset at 325ms before preflashevent

% preprocess LFPs
Fs = D.lfpFs;
D.adjLfpsClean = D.adjLfps; %interpolateLfpOverSpikeTimes(D.adjLfps, channelsToLoad, Fs, D.allMUAStructs);
hiCutoffFreq = 100;
[channelDataCARNorm,channelDataNorm,commonAverageNorm,isNoisyChannel] = preprocessLfps(...
        D.adjLfpsClean, Fs, D.lfpNames, processedDataDir, [], hiCutoffFreq, 1, []);
D.adjLfpsClean = [];
flashEventsClean = detectOutlierLfpEvents(channelDataCARNorm, Fs, D.lfpNames, origFlashEvents, ...
        postFlashWindowOffset, processedDataDir, [], 1, []);
preFlashEventsCleanWindowOffset = [min([baselineWindowOffset saccadeWindowOffset]) max([baselineWindowOffset saccadeWindowOffset])];
preFlashEventsClean = detectOutlierLfpEvents(channelDataCARNorm, Fs, D.lfpNames, preFlashesEvents, ...
        preFlashEventsCleanWindowOffset, processedDataDir, [], 1, []);
   
nFlashes = numel(flashEventsClean);
nTrials = numel(preFlashEventsClean);
startIndicesStim = round((flashEventsClean + postFlashWindowOffset(1)) * Fs); % time to index conversion
endIndicesStim = startIndicesStim + round(diff(postFlashWindowOffset) * Fs) - 1;
startIndicesBl = round((preFlashEventsClean + baselineWindowOffset(1)) * Fs); % time to index conversion
endIndicesBl = startIndicesBl + round(diff(baselineWindowOffset) * Fs) - 1;
startIndicesSac = round((preFlashEventsClean + saccadeWindowOffset(1)) * Fs); % time to index conversion
endIndicesSac = startIndicesSac + round(diff(saccadeWindowOffset) * Fs) - 1;
t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
nTime = numel(t);
% assert(all(nTime == (endIndicesStim - startIndicesStim + 1)));

nChannels = length(channelsToLoad);
responses = nan(nChannels, nTime, nFlashes);
baseline = nan(nChannels, nTime, nTrials);
saccade = nan(nChannels, nTime, nTrials);
for j = 1:nChannels
    % can vectorize???
    for i = 1:nFlashes
        responses(j,:,i) = channelDataCARNorm(j,startIndicesStim(i):endIndicesStim(i));
    end
    for k = 1:nTrials
        baseline(j,:,k) = channelDataCARNorm(j,startIndicesBl(k):endIndicesBl(k));
        saccade(j,:,k) = channelDataCARNorm(j,startIndicesSac(k):endIndicesSac(k));
    end
end
allResponses = cat(3,responses,saccade);

% create covariance matrices
responsesCont = reshape(allResponses, nChannels, []);
responsesCont = bsxfun(@minus, responsesCont, mean(responsesCont,2));
responsesCov = (responsesCont*responsesCont') / size(responsesCont,2);
baselineCont = reshape(baseline, nChannels, []);
baselineCont = bsxfun(@minus, baselineCont, mean(baselineCont,2));
baselineCov = (baselineCont*baselineCont') / size(baselineCont,2);
[evecs, evals] = eig(responsesCov,baselineCov);
[~,eigidx] = sort(diag(evals));
evecs = evecs(:,eigidx);

responsesGed = reshape( (responsesCont' * evecs)',nChannels,nTime,nFlashes+nTrials);
responsesGedCompMaps = reshape( (responsesCov' * evecs(:,1))',1,nChannels,1);
responsesGedCompMaps2 = reshape( (responsesCov' * evecs(:,2))',1,nChannels,1);

figure;
plot(responsesGedCompMaps,1:32)
set(gca, 'YDir', 'reverse');

figure;
plot(responsesGedCompMaps2,1:32)
set(gca, 'YDir', 'reverse');

baselineGed = reshape( (baselineCont' * evecs)',nChannels,nTime,nTrials);

% plot it
averageResponse = mean(allResponses, 3); % average across flashes
averageResponseGED = mean(responsesGed, 3); % average across flashes
responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(averageResponse)));
figure;
subplot(221)
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset);
%     plot(t, responsePlotYScale*averageResponseGED(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset) - 1]);
    end
end
subplot(222)
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
%     plot(t, responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset);
    plot(t, responsePlotYScale*averageResponseGED(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*averageResponseGED(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*averageResponseGED(j,:) - j*responsePlotYOffset) - 1]);
    end
end

baselineResponse = mean(baseline, 3); % average across flashes
baselineResponseGED = mean(baselineGed, 3); % average across flashes
subplot(223)
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*baselineResponse(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*baselineResponse(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*baselineResponse(j,:) - j*responsePlotYOffset) - 1]);
    end
end
subplot(224)
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*baselineResponseGED(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*baselineResponseGED(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*baselineResponseGED(j,:) - j*responsePlotYOffset) - 1]);
    end
end

[tBaseline,averageBaselineResponses] = computeResponsesInWindow(channelDataCARNorm, preFlashesEvents, ...
        baselineWindowOffset, Fs);
[tBaselineExp,averageBaselineResponsesExpCAR] = computeResponsesInWindow(channelDataCARNorm, preFlashesEvents, ...
        expandedPlotWindowOffset, Fs);

    
%% save responses to mat file
saveFileName = sprintf('%s/%s-%s-responses-v%d.mat', processedDataDir, plotFileNamePrefix, ref, v);
save(saveFileName, 'D', 'R', 'responses', 't', 'periFlashWindowOffset', 'isNoisyChannel');

%% subtract out baseline
% units are standard deviations from baselined mean
averageResponse = mean(responses, 3); % average across flashes

% subtract mean baseline activity (-0.1, 0] seconds before flash
flashBaselineWindowOffset = [-0.1 0];

% TODO make so that baseline can have arbitrary end time
indexFlashTime = -round(postFlashWindowOffset(1) * Fs);
indexStartBaselineFlashTime = -round((postFlashWindowOffset(1) - flashBaselineWindowOffset(1)) * Fs) + 1;
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
xlim(postFlashWindowOffset);
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from Flash Onset (s)');
title(sprintf('%s %s - Response to Full-Field Flash (N=%d) (%s)', sessionName, areaName, nFlashes, ref));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, minYLim+1, '35-45 ms', 'Color', 'm');

plotFileName = sprintf('%s/%s-%s-lfpLines-v%d.png', processedDataDir, plotFileNamePrefix, ref, v);
%export_fig(plotFileName, '-nocrop');

%% staggered channel color plot of average visually evoked LFP
figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.1);
hold on;
imagesc(t, 1:nChannels, averageResponse);
set(gca, 'YDir', 'reverse');
plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim([0.5 nChannels+0.5]);
xlabel('Time from Flash Onset (s)');
ylabel('Channel Number (1 = topmost)');
title(sprintf('%s %s - Response to Full-Field Flash (N=%d) (%s)', sessionName, areaName, nFlashes, ref));
maxCAxis = max(abs(caxis));
caxis([-maxCAxis maxCAxis]);
colormap(getCoolWarmMap());
colorbar;
set(gca, 'FontSize', 16);

if (nChannels > 31 && strcmp(ref, 'BIP')) || (nChannels > 32 && strcmp(ref, 'CAR'))
    plot(xlim(), [nChannels/2+0.5 nChannels/2+0.5], 'Color', 0.3*ones(3, 1));
end

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, nChannels, '35-45 ms', 'Color', 'm');

plotFileName = sprintf('%s/%s-%s-lfpColor-v%d.png', processedDataDir, plotFileNamePrefix, ref, v);
%export_fig(plotFileName, '-nocrop');

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

plotFileName = sprintf('%s/%s-%s-lfpColorBounds-v%d.png', processedDataDir, plotFileNamePrefix, ref, v);
%export_fig(plotFileName, '-nocrop');

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
xlabel(sprintf('Mean SDs Early Response to Flash (%d-%d ms)', earlyActivityWindowOffset * 1000));
ylabel('Channel Number');
grid on;
ylim(channelIndices([1 end]) + [-1 1]);
xlim([-1 1] * max(abs(xlim())));
title(sprintf('%s %s - Mean Early Response (%s)', sessionName, areaName, ref));

plotFileName = sprintf('%s/%s-%s-meanEarlyResponse-v%d.png', processedDataDir, plotFileNamePrefix, ref, v);
%export_fig(plotFileName, '-nocrop');





%% staggered channel color plot of average visually evoked LFP
figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.1);
hold on;
imagesc(t, 1:nChannels, averageResponse);
set(gca, 'YDir', 'reverse');
plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim([0.5 nChannels+0.5]);
xlabel('Time from Flash Onset (s)');
ylabel('Channel Number (1 = topmost)');
title(sprintf('%s %s - Response to Full-Field Flash (N=%d) (%s)', sessionName, areaName, nFlashes, ref));
maxCAxis = max(abs(caxis));
caxis([-maxCAxis maxCAxis]);
colormap(getCoolWarmMap());
colorbar;
set(gca, 'FontSize', 16);

if (nChannels > 31 && strcmp(ref, 'BIP')) || (nChannels > 32 && strcmp(ref, 'CAR'))
    plot(xlim(), [nChannels/2+0.5 nChannels/2+0.5], 'Color', 0.3*ones(3, 1));
end

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, nChannels, '35-45 ms', 'Color', 'm');

plotFileName = sprintf('%s/%s-%s-lfpColor-v%d.png', processedDataDir, plotFileNamePrefix, ref, v);
%export_fig(plotFileName, '-nocrop');





% eof