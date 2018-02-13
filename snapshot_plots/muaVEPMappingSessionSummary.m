
%% write updated D struct
% remove wf and ts to save space
for j = 1:nUnits
    D.allMUAStructs{j}.wf = [];
    D.allMUAStructs{j}.ts = [];
end
saveFileName = sprintf('%s/D-afterSpikeVEPM-%s.mat', processedDataDir, blockName);
fprintf('Writing post-processed data file: %s ...', saveFileName);
tic;
save(saveFileName, 'D');
fprintf(' done (%0.2f s).\n', toc);

%% plot cdf of latencies
nUnits = numel(D.allMUAStructs);
fprintf('Processing %d units...\n', nUnits);

latencies = nan(nUnits, 3); % latency, channel number, spike struct index
for j = 1:nUnits
    spikeStruct = D.allMUAStructs{j};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
            nUnits, round(j/nUnits*100));
    
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

plotFileName = sprintf('%s/%s_%s-latencies-%s.png', processedDataDir, sessionName, areaName, blockName);
export_fig(plotFileName, '-nocrop');

%% scatter plot of latencies by channel

% get mean latency by channel
% assume channel numbering starts at 1
meanLatenciesByChannel = nan(D.nMUACh, 1);
for i = 1:numel(meanLatenciesByChannel)
    meanLatenciesByChannel(i) = nanmean(latencies(latencies(:,2) == R.muaChannelsToLoad(i)), 1);
end

channelIndices = R.muaChannelsToLoad;
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

plotFileName = sprintf('%s/%s_%s-latenciesByChannel-%s.png', processedDataDir, sessionName, areaName, blockName);
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
for j = 1:nUnits
    spikeStruct = D.allMUAStructs{j};
    unitName = spikeStruct.name;
    nameParts = strfind(unitName, '_');
    unitNameShort = unitName((nameParts(end)+1):end);
    spikeTimes = spikeStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
            nUnits, round(j/nUnits*100));
    
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

plotFileName = sprintf('%s/%s_%s-spdfByChannelJoyplotFilled-all-%s.png', processedDataDir, sessionName, areaName, blockName);
export_fig(plotFileName, '-nocrop');

%% joyplot of spdfs by channel, scaled by SDs of baseline, cells only

figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.10, 'MB', 0.10, 'MR', 0.05);
hold on;
plotCount = 0;
yScale = 1/15;
col = lines(1);
cmap = winter();
cmap = flipud(cmap(floor(size(cmap, 1)/2):end,:));
% cmap = getCoolWarmMap();
for j = 1:nUnits
    spikeStruct = D.allMUAStructs{j};
    unitName = spikeStruct.name;
    nameParts = strfind(unitName, '_');
    unitNameShort = unitName((nameParts(end)+1):end);
    spikeTimes = spikeStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
            nUnits, round(j/nUnits*100));
    
    % include all units even if there was insignificant response
    % but exclude if no firing or super sparse firing
    if any(isnan(spikeStruct.vepmPsthParams.normPsthResponse)) || ...
            spikeStruct.vepmPsthParams.meanResponseByPsth < minResponseByPsth
        continue;
    end
    
    % exclude if axon
    if ~(strcmp(spikeStruct.physClass, 'Broad-Spiking') || strcmp(spikeStruct.physClass, 'Narrow-Spiking'))
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

plotFileName = sprintf('%s/%s_%s-spdfByChannelJoyplotFilled-cellsOnly-%s.png', processedDataDir, sessionName, areaName, blockName);
export_fig(plotFileName, '-nocrop');

%% joyplot of spdfs by channel, scaled by SDs of baseline, cells only

figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.10, 'MB', 0.10, 'MR', 0.05);
hold on;
plotCount = 0;
yScale = 1/10;
cols = lines(2);
bsCol = cols(1,:);
nsCol = cols(2,:);
cmap = gray();
cmap = flipud(cmap(floor(size(cmap, 1)/2):end,:));
% cmap = getCoolWarmMap();
for j = 1:nUnits
    spikeStruct = D.allMUAStructs{j};
    unitName = spikeStruct.name;
    nameParts = strfind(unitName, '_');
    unitNameShort = unitName((nameParts(end)+1):end);
    spikeTimes = spikeStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
            nUnits, round(j/nUnits*100));
    
    % include all units even if there was insignificant response
    % but exclude if no firing or super sparse firing
    if any(isnan(spikeStruct.vepmPsthParams.normPsthResponse)) || ...
            spikeStruct.vepmPsthParams.meanResponseByPsth < minResponseByPsth
        continue;
    end
    
    % exclude if axon
    if ~(strcmp(spikeStruct.physClass, 'Broad-Spiking') || strcmp(spikeStruct.physClass, 'Narrow-Spiking'))
        continue;
    end
    
    if strcmp(spikeStruct.physClass, 'Broad-Spiking')
        col = bsCol;
    else
        col = nsCol;
    end
    
    plotCount = plotCount + 1;
    % reverse yVal since we reverse the axis
    yVal = -1 * spikeStruct.vepmPsthParams.normPsthResponse * yScale + plotCount;
    plot(spikeStruct.vepmPsthParams.t, yVal, '-', 'Color', col, 'LineWidth', 1);
%     fillX = [spikeStruct.vepmPsthParams.t spikeStruct.vepmPsthParams.t([end 1])];
%     fillY = [yVal plotCount plotCount];
%     fillC = [spikeStruct.vepmPsthParams.normPsthResponse 0 0];
%     fill(fillX, fillY, fillC, 'FaceAlpha', 0.85);
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
title(sprintf('%s %s - Flash Responses, Cells Only (N=%d)', sessionName, areaName, plotCount));
legendTextBS = sprintf('{\\color[rgb]{%f,%f,%f}%s}', bsCol, 'Broad-Spiking');
legendTextNS = sprintf('{\\color[rgb]{%f,%f,%f}%s}', nsCol, 'Narrow-Spiking');
text(0.85, -0.07, {legendTextBS, legendTextNS}, 'FontSize', 10, 'Units', 'normalized');

plotFileName = sprintf('%s/%s_%s-spdfByChannelJoyplot_cellsOnly-%s.png', processedDataDir, sessionName, areaName, blockName);
export_fig(plotFileName, '-nocrop');

%% joyplot of spdfs by channel, scaled by SDs of baseline, cells only

figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.10, 'MB', 0.10, 'MR', 0.05);
hold on;
plotCount = 0;
yScale = 1/10;
cols = lines(5);
bsCol = cols(1,:);
nsCol = cols(2,:);
otherCol = cols(5,:);
cmap = gray();
cmap = flipud(cmap(floor(size(cmap, 1)/2):end,:));
% cmap = getCoolWarmMap();
for j = 1:nUnits
    spikeStruct = D.allMUAStructs{j};
    unitName = spikeStruct.name;
    nameParts = strfind(unitName, '_');
    unitNameShort = unitName((nameParts(end)+1):end);
    spikeTimes = spikeStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
            nUnits, round(j/nUnits*100));
    
    % include all units even if there was insignificant response
    % but exclude if no firing or super sparse firing
    if any(isnan(spikeStruct.vepmPsthParams.normPsthResponse)) || ...
            spikeStruct.vepmPsthParams.meanResponseByPsth < minResponseByPsth
        continue;
    end
    
    if strcmp(spikeStruct.physClass, 'Broad-Spiking')
        col = bsCol;
    elseif strcmp(spikeStruct.physClass, 'Narrow-Spiking')
        col = nsCol;
    else
        col = otherCol;
    end
    
    plotCount = plotCount + 1;
    % reverse yVal since we reverse the axis
    yVal = -1 * spikeStruct.vepmPsthParams.normPsthResponse * yScale + plotCount;
    plot(spikeStruct.vepmPsthParams.t, yVal, '-', 'Color', col, 'LineWidth', 1);
%     fillX = [spikeStruct.vepmPsthParams.t spikeStruct.vepmPsthParams.t([end 1])];
%     fillY = [yVal plotCount plotCount];
%     fillC = [spikeStruct.vepmPsthParams.normPsthResponse 0 0];
%     fill(fillX, fillY, fillC, 'FaceAlpha', 0.85);
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
legendTextBS = sprintf('{\\color[rgb]{%f,%f,%f}%s}', bsCol, 'Broad-Spiking');
legendTextNS = sprintf('{\\color[rgb]{%f,%f,%f}%s}', nsCol, 'Narrow-Spiking');
legendTextOther = sprintf('{\\color[rgb]{%f,%f,%f}%s}', otherCol, 'Other');
text(0.85, -0.07, {legendTextBS, legendTextNS, legendTextOther}, 'FontSize', 10, 'Units', 'normalized');

plotFileName = sprintf('%s/%s_%s-spdfByChannelJoyplot_allColorCoded-%s.png', processedDataDir, sessionName, areaName, blockName);
export_fig(plotFileName, '-nocrop');

