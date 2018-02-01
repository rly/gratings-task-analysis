
clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';

sessionInd = 8;

%% load recording information
recordingInfo = readRecordingInfo();
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% setup and load data
fprintf('\n-------------------------------------------------------\n');
fprintf('Unit Analysis\n');
fprintf('Loading %s...\n', pl2FilePath);

tic;
isLoadSpikes = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
isLoadDirect = 0;
D = loadPL2(pl2FilePath, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad); 

processedDataDir = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

%% waveform analysis
nUnits = numel(D.allSpikeStructs);
fprintf('Processing %d units...\n', nUnits);
assert(numel(blockNames) == numel(D.blockStartTimes));

%% plot waveforms for each cell over time
yScale = 50;
col = [0 0 0];
for j = 1:nUnits
    spikeStruct = D.allSpikeStructs{j};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
            nUnits, round(j/nUnits*100));

    nTimeWaveform = size(D.allSpikeStructs{1}.wf, 2);
    meanSpikeWaveformByBlock = nan(numel(blockNames), nTimeWaveform);
    seSpikeWaveformByBlock = nan(numel(blockNames), nTimeWaveform);
    spikeCountsByBlock = nan(numel(blockNames), 1);
    spikeRatesByBlock = nan(numel(blockNames), 1);
    for k = 1:numel(blockNames)
        blockStartTime = D.blockStartTimes(k);
        blockStopTime = D.blockStopTimes(k);
        
        spikeIndices = D.allSpikeStructs{j}.ts >= blockStartTime & ...
                D.allSpikeStructs{j}.ts <= blockStopTime;

        wfs = D.allSpikeStructs{j}.wf(spikeIndices,:);
        if ~isempty(wfs)
            meanSpikeWaveformByBlock(k,:) = mean(wfs);
            seSpikeWaveformByBlock(k,:) = std(wfs) / sqrt(size(wfs, 1));
        end
        spikeCountsByBlock(k) = sum(spikeIndices);
        spikeRatesByBlock(k) = spikeCountsByBlock(k) / (blockStopTime - blockStartTime);
    end
    
    t = ((0:size(meanSpikeWaveformByBlock, 2)-1) / D.timestampFrequency - D.allSpikeStructs{1}.thresholdTime) * 1000;
    
    %% test for significant spike rate drift over the gratings task blocks
    firstGTaskBlockStartTime = D.blockStartTimes(gratingsTask3DIndices(1));
    lastGTaskBlockStopTime = D.blockStopTimes(gratingsTask3DIndices(end));
    windowSize = 300; % seconds
    movingWindowStep = floor(windowSize/2);
    
    spikeIndices = D.allSpikeStructs{j}.ts >= firstGTaskBlockStartTime & ...
                D.allSpikeStructs{j}.ts <= lastGTaskBlockStopTime;
    spikeTimes = D.allSpikeStructs{j}.ts(spikeIndices);
    
    binStartTimes = firstGTaskBlockStartTime:movingWindowStep:lastGTaskBlockStopTime-windowSize;
    numWindows = numel(binStartTimes);
    spikesPerBin = nan(numWindows, 1);
    for k = 1:numWindows
        spikesPerBin(k) = sum(spikeTimes >= binStartTimes(k) & spikeTimes < binStartTimes(k) + windowSize);
    end
    
    % test for time series stationarity
    % one issue: this method preferentially labels decreasing time series
    % as stationary. reversing the order of the data yields different
    % results.
    [hIsSpikeRateStationaryFwd,pIsSpikeRateStationaryFwd,stat,cValue,reg] = adftest(spikesPerBin, 'model', 'TS', 'alpha', 0.1);
    [hIsSpikeRateStationaryRev,pIsSpikeRateStationaryRev,stat,cValue,reg] = adftest(spikesPerBin(end:-1:1), 'model', 'TS', 'alpha', 0.1);
    
    % test for monotonic trend with strong alpha
    % REALLY SHOULD be using the reverse test where H0 = there is a trend
    % or the above stationarity test, perhaps with a bigger window size
    [hIsSpikeRateTrend,pIsSpikeRateTrend] = Mann_Kendall(spikesPerBin, 0.005);
    
    minMeanSpikeRate = 0.2;
    if ~hIsSpikeRateTrend && mean(spikesPerBin) / windowSize >= minMeanSpikeRate
        fprintf('Hooray! Null hypothesis not rejected...\n');
    end
    
    if (hIsSpikeRateStationaryFwd && hIsSpikeRateStationaryRev) && mean(spikesPerBin) / windowSize >= minMeanSpikeRate
        fprintf('Hooray! Spike rate is stationary.\n');
    end
    
    % alternatively, test whether bins are drawn from a uniform
    % distribution
    
    figure;
    plot(spikesPerBin);
%     if hIsSpikeRateStationary
%         title('Hooray! Spike rate is stationary.');
%         fprintf('Hooray! Spike rate is stationary.\n');
%         stop
%     else
%         title('Spike rate has a unit root. It is not stationary.');
%     end
    
    
        
    
    %% create figure
    f = figure_tr_inch(13, 7.5); clf;
    set(gcf, 'Color', 'white');
    set(gcf, 'renderer', 'painters');

    %% make main title
    axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
    set(get(axBig, 'Title'), 'Visible', 'on')

    modTitle = sprintf('Spike Analysis by Block: %s', unitName);
    titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
    title(modTitle, 'FontSize', 14, titleParams{:});

    %% location params
    waveformW = 0.1;
    bigWaveformW = 0.36;
    rateHistW = 0.36;
    waveformH = 0.15;
    bigWaveformH = 0.84;
    rateHistH = bigWaveformH;

    waveformLeft = 0.05;
    bigWaveformLeft = waveformLeft + waveformW + 0.08;
    rateHistLeft = bigWaveformLeft + bigWaveformW + 0.03;

    btm = 0.07;
    infoTextTop = btm + 0.1;
    infoText2Top = infoTextTop + 0.53;
    waveformBtm = 0.76;
    
    %% plot spike waveform
    axes('Position', [waveformLeft waveformBtm waveformW waveformH]); 
    plotSpikeWaveform(D, j);

    %% info
    writeUnitInfo(spikeStruct, axBig, -0.03, infoTextTop);
    
    %% waveform across blocks
    axes('Position', [bigWaveformLeft btm bigWaveformW bigWaveformH]); 
    hold on;
    for k = 1:numel(blockNames)
        if ~all(isnan(meanSpikeWaveformByBlock(k,:)))
            yVal = meanSpikeWaveformByBlock(k,:) * yScale + k;
            jbfill(t, yVal + seSpikeWaveformByBlock(k,:) * yScale, ...
                    yVal - seSpikeWaveformByBlock(k,:) * yScale, col, col, 0.3);
            hold on;
            plot(t, yVal, '-', 'Color', col, 'LineWidth', 2);
        else
            plot(t, zeros(size(t)) + k, '--', 'Color', col);
        end
    end
    xlim([t(1) t(end)]);
    ylim([-2 numel(blockNames) + 3]);
    xlabel('Time from Threshold (ms)');
    ylabel('Block Name');
    title('Waveform', 'Interpreter', 'none');
    grid on;
    box off;
    set(gca, 'XTick', -0.4:0.2:1);
    set(gca, 'YTick', 1:numel(blockNames));
    set(gca, 'YTickLabel', blockNames);
    text(0.98, 0.02, 'Mean +/- 1 SEM', 'FontSize', 8, 'HorizontalAlignment', 'right', 'Units', 'normalized');

    %% firing rate by block
    axes('Position', [rateHistLeft btm rateHistW rateHistH]); 
    hold on;
    barCol = 0.3 * ones(3, 1);
    barh(spikeRatesByBlock, 'FaceColor', barCol, 'EdgeColor', barCol);
    ylim([-2 numel(blockNames) + 3]);
    xlabel('Firing Rate (Hz)');
    title('Firing Rate', 'Interpreter', 'none');
    set(gca, 'YTickLabel', []);
    set(gca, 'YTick', []);
    
    %% save
    plotFileName = sprintf('%s/%s-all-unitAnalysisByBlock.png', processedDataDir, unitName);
    export_fig(plotFileName, '-nocrop');
    close;
    
    %% make bar plot where time is on x axis, spike rates are plotted every 2 minutes
    % and shaded areas show when there is a block, with text asaying the
    % block name
end

%% check refractory period violations (spikes within 2 ms of each other)

%% histogram of trough to peak times
isCell = false(nUnits, 1);
isBroadSpiking = false(nUnits, 1);
isNarrowSpiking = false(nUnits, 1);
troughToPeakTimes = nan(nUnits, 1);

for j = 1:nUnits
    spikeStruct = D.allSpikeStructs{j};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    
    isCell(j) = strcmp(spikeStruct.physClass, 'Broad-Spiking') || strcmp(spikeStruct.physClass, 'Narrow-Spiking');
    if isCell(j)
        isBroadSpiking(j) = strcmp(spikeStruct.physClass, 'Broad-Spiking');
        isNarrowSpiking(j) = strcmp(spikeStruct.physClass, 'Narrow-Spiking');
    end
    
    % could require a minimum firing rate but this is across the whole task
    % and cells fade in and out, warping the firing rate
    % could do a rolling window and use the maximum firing rate across all
    % windows
    
    troughToPeakTimes(j) = spikeStruct.troughToPeakTime * 1000;
end

figure_tr_inch(8, 6); clf;
subaxis(1, 1, 1, 'ML', 0.08, 'MR', 0.03);
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');
hold on;
x = 0:0.02499:0.8; % if we use 0.025 we run into rounding issues
histogram(troughToPeakTimes(isBroadSpiking), x);
histogram(troughToPeakTimes(isNarrowSpiking), x);
ylabel('Number of Cells');
xlabel('Trough to Peak Time (ms)');
box off;
grid on;
set(gca, 'FontSize', 12);
title(sprintf('%s - %s (N=%d)', sessionName, areaName, sum(isCell)), 'FontSize', 16);

plotFileName = sprintf('%s/all-troughToPeakTimeHist.png', processedDataDir);
export_fig(plotFileName, '-nocrop');
% close;

%% plot BS vs NS cells by channel/depth



