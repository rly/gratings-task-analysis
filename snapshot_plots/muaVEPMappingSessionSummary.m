function muaVEPMappingSessionSummary(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad)
% Session Summary of MUA VEP Mapping

%% setup and load data
v = 11;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task MUA Analysis\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('MUA Channels to Load: %d\n', muaChannelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

nUnits = numel(muaChannelsToLoad);
muaChannelRangeStr = sprintf('%d_%d', muaChannelsToLoad([1 end]));

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);
recordingInfo = rmfield(recordingInfo, 'muaChannelsToLoad'); % remove b/c we're using passed value
R = recordingInfo(sessionInd);
sessionName = R.sessionName;
areaName = R.areaName;

scriptName = 'MUA_VEPM';
processedDataDir = sprintf('%s/%s/%s', processedDataRootDir, sessionName, scriptName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end

blockIndices = R.vepmIndices;
blockName = strjoin(R.blockNames(blockIndices), '-');

%% load data from individual mat files
allMUAStructs = cell(nUnits, 1);
for j = 1:nUnits
    unitName = sprintf('%s_%s_%dM', sessionName, areaName, muaChannelsToLoad(j)); % ideally find a better way
    saveFileName = sprintf('%s/%s-%s-vepm-v%d.mat', processedDataDir, unitName, blockName, v);
    fprintf('(%d/%d = %d%%) Loading %s...\n', j, nUnits, round(j/nUnits * 100), saveFileName);
    S = load(saveFileName, 'muaStruct');
    allMUAStructs{j} = S.muaStruct;
end

latencyWindowOffset = allMUAStructs{1}.vepmPsthParams.latencyWindowOffset;

%% plot cdf of latencies
fprintf('Processing %d units...\n', nUnits);

latencies = nan(nUnits, 3); % latency, channel number, spike struct index
for j = 1:nUnits
    spikeStruct = allMUAStructs{j};
    
    if isstruct(spikeStruct.latencyInfo)
        latencies(j,1) = spikeStruct.latencyInfo.latency * 1000;
        latencies(j,2) = spikeStruct.channelID;
        latencies(j,3) = j;
    end
end
latencies = trimNanRows(latencies);

figure_tr_inch(6, 6);
subaxis(1, 1, 1, 'ML', 0.15, 'MB', 0.13, 'MR', 0.07);
plotH = cdfplot(latencies(:,1));
set(plotH, 'LineWidth', 3);
set(gca, 'FontSize', 16);
set(gca, 'XTick', 25:25:225);
set(gca, 'YTick', 0:0.25:1);
set(gca, 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1);
xlim(latencyWindowOffset * 1000);
xlabel('Response Latency to Flash (ms)');
ylabel('Proportion of Units');
title(sprintf('%s %s - Flash Latencies (N=%d)', sessionName, areaName, size(latencies, 1)));

plotFileName = sprintf('%s/%s_%s_%s-latencies-%s-v%d.png', processedDataDir, sessionName, areaName, muaChannelRangeStr, blockName, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% scatter plot of latencies by channel

% get mean latency by channel
% assume channel numbering starts at 1
meanLatenciesByChannel = nan(numel(muaChannelsToLoad), 1);
for i = 1:numel(meanLatenciesByChannel)
    meanLatenciesByChannel(i) = nanmean(latencies(latencies(:,2) == muaChannelsToLoad(i)), 1);
end

channelIndices = muaChannelsToLoad;
channelIndices(isnan(meanLatenciesByChannel)) = [];
meanLatenciesByChannel(isnan(meanLatenciesByChannel)) = [];

figure_tr_inch(6, 6);
subaxis(1, 1, 1, 'ML', 0.15, 'MB', 0.13, 'MR', 0.07);
hold on;
plot(latencies(:,1), latencies(:,2), '.', 'MarkerSize', 25);
plot(meanLatenciesByChannel, channelIndices, '--', 'LineWidth', 2, 'Color', lines(1));
set(gca, 'FontSize', 16);
set(gca, 'YDir', 'reverse');
xlabel('Response Latency to Flash (ms)');
ylabel('Channel Number');
grid on;
ylim(channelIndices([1 end]) + [-1 1]);
title(sprintf('%s %s - Flash Latencies (N=%d)', sessionName, areaName, numel(meanLatenciesByChannel)));

plotFileName = sprintf('%s/%s_%s_%s-latenciesByChannel-%s-v%d.png', processedDataDir, sessionName, areaName, muaChannelRangeStr, blockName, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% joyplot of spdfs by channel, scaled by SDs of baseline

minResponseByPsth = 0.2; % Hz

figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.10, 'MB', 0.10, 'MR', 0.05);
hold on;
plotCount = 0;
yScale = 1/15;
col = lines(1);
cmap = winter();
cmap = flipud(cmap(floor(size(cmap, 1)/2):end,:));
% cmap = getCoolWarmMap();
xBounds = [-0.2 0.27];
for j = 1:nUnits
    spikeStruct = allMUAStructs{j};
    unitName = spikeStruct.name;
    nameParts = strfind(unitName, '_');
    unitNameShort = unitName((nameParts(end)+1):end);
    
    % include all units even if there was insignificant response
    % but exclude if no firing or super sparse firing
    if any(isnan(spikeStruct.vepmPsthParams.normPsthResponse)) || ...
            spikeStruct.vepmPsthParams.meanResponseByPsth < minResponseByPsth
        continue;
    end
    
    plotCount = plotCount + 1;
    % reverse yVal since we reverse the axis
    yVal = -1 * spikeStruct.vepmPsthParams.normPsthResponse * yScale + plotCount;
    plot(spikeStruct.vepmPsthParams.t, yVal, '-', 'Color', col);
    fillX = [spikeStruct.vepmPsthParams.t spikeStruct.vepmPsthParams.t([end 1])];
    fillY = [yVal plotCount plotCount];
    fillC = [spikeStruct.vepmPsthParams.normPsthResponse 0 0];
    fill(fillX, fillY, fillC, 'FaceAlpha', 0.85);
    text(xBounds(1) - range(xBounds)/100, plotCount, unitNameShort, 'FontSize', 10, 'HorizontalAlignment', 'right');
end
cBounds = max(abs(caxis)) * [-1 1];
caxis(cBounds); % symmetric cmap
colormap(cmap);
plot([0 0], [0 plotCount+1], '-', 'Color', 0.3*ones(3, 1)); 
set(gca, 'FontSize', 16);
set(gca, 'YDir', 'reverse');
xlabel('Time from Flash Onset (s)');
ylabel({'Unit (Ordered by Depth)', ''});
set(gca, 'YTickLabel', '');
grid on;
xlim(xBounds);
ylim([0 plotCount + 1]);
title(sprintf('%s %s - Flash Responses (N=%d)', sessionName, areaName, plotCount));

plotFileName = sprintf('%s/%s_%s_%s-spdfByChannelJoyplotFilled-all-%s-v%d.png', ...
        processedDataDir, sessionName, areaName, muaChannelRangeStr, blockName, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% heatmap
figure_tr_inch(4, 5.5);
ax1 = subaxis(1, 1, 1, 'ML', 0.07, 'MB', 0.14, 'MR', 0.06, 'MT', 0.03);
hold on;
cmap = getCoolWarmMap();
xBounds = [-0.05 0.251];

spikeStruct = allMUAStructs{1};
t = spikeStruct.vepmPsthParams.t;
normPsthResp = nan(nUnits, numel(spikeStruct.vepmPsthParams.normPsthResponse));
for j = 1:nUnits
    spikeStruct = allMUAStructs{j};
    normPsthResp(j,:) = spikeStruct.vepmPsthParams.normPsthResponse;
end
imagesc(t, 1:nUnits, normPsthResp);

cBounds = max(abs(caxis)) * [-1 1] * 1/2;
caxis(cBounds); % symmetric cmap
colormap(cmap);
plot([0 0], [0.5 nUnits+0.5], '-', 'Color', 0.3*ones(3, 1)); 
set(gca, 'FontSize', 16);
set(gca, 'YDir', 'reverse');
xlabel('Time from Flash Onset (s)');
% ylabel({'Unit (Ordered by Depth)', ''});
% set(gca, 'XTickLabel', {' '});
set(gca, 'XTick', -0.05:0.05:0.25);
set(gca, 'YTickLabel', '');
grid on;
xlim(xBounds);
ylim([0.5 nUnits+0.5]);
ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none', 'FontSize', 16);
set(ax2, 'XLim', get(ax1, 'XLim'), 'YLim', get(ax1, 'YLim'));
set(ax2, 'YTick', [], 'YColor', 'w', 'YAxisLocation', 'right');
set(ax2, 'XTick', get(ax1, 'XTick'), 'XTickLabel', [], 'XAxisLocation', 'bottom', 'TickDir', 'out');

% title(sprintf('%s %s - Flash Responses (N=%d)', sessionName, areaName, nUnits));

plotFileName = sprintf('%s/%s_%s_%s-spdfByChannelHeat-all-%s-v%d.png', ...
        processedDataDir, sessionName, areaName, muaChannelRangeStr, blockName, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');


