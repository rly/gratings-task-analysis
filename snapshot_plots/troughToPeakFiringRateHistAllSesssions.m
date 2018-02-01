% across sessions

clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';

nUnitsApprox = 100; % make sure this is an underestimate
isCell = false(nUnitsApprox, 1);
isBroadSpiking = false(nUnitsApprox, 1);
isNarrowSpiking = false(nUnitsApprox, 1);
troughToPeakTimes = nan(nUnitsApprox, 1);
meanFiringRates = nan(nUnitsApprox, 1);
meanWfs = nan(nUnitsApprox, 56);
localization = cell(nUnitsApprox, 1);
unitCount = 0;
minFiringRate = 0.5;

fprintf('\n-------------------------------------------------------\n');
fprintf('Waveform Trough to Peak Time and Firing Rate Analysis\n');

sessionIndAll = [1 2 3 4 5 6 8 9];
for k = 1:numel(sessionIndAll)
    sessionInd = sessionIndAll(k);

    %% load recording information
    recordingInfo = readRecordingInfo();
    struct2var(recordingInfo(sessionInd));
    pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

    %% setup and load data
    fprintf('Loading %s...\n', pl2FilePath);

    tic;
    isLoadSpikes = 1;
    isLoadLfp = 0;
    isLoadSpkc = 0;
    isLoadDirect = 0;
    D = loadPL2(pl2FilePath, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc, isLoadDirect, ...
            spikeChannelPrefix, spikeChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad); 

    fprintf('... done (%0.2f s).\n', toc);

    nUnits = numel(D.allSpikeStructs);
    fprintf('Processing %d units...\n', nUnits);
    
    %% histogram of trough to peak times
    for j = 1:nUnits
        spikeStruct = D.allSpikeStructs{j};
        unitName = spikeStruct.name;
        spikeTimes = spikeStruct.ts;
        unitCount = unitCount + 1;
        
        if ismember(spikeStruct.channelID, pldChannels)
            localization{unitCount} = 'PLd';
        elseif ismember(spikeStruct.channelID, plvChannels)
            localization{unitCount} = 'PLv';
        elseif ismember(spikeStruct.channelID, pmChannels)
            localization{unitCount} = 'PM';
        elseif ismember(spikeStruct.channelID, piChannels)
            localization{unitCount} = 'PI';
        else
            localization{unitCount} = '';
        end

        isCell(unitCount) = strcmp(spikeStruct.physClass, 'Broad-Spiking') || strcmp(spikeStruct.physClass, 'Narrow-Spiking');
        if isCell(unitCount)
            isBroadSpiking(unitCount) = strcmp(spikeStruct.physClass, 'Broad-Spiking');
            isNarrowSpiking(unitCount) = strcmp(spikeStruct.physClass, 'Narrow-Spiking');
        else
            isBroadSpiking(unitCount) = false;
            isNarrowSpiking(unitCount) = false;
        end

        troughToPeakTimes(unitCount) = spikeStruct.troughToPeakTime * 1000;
        meanWfs(unitCount,:) = spikeStruct.meanWf;
        
        % calculate the mean firing rate between the 25th and 75th
        % percentile of spike times, regardless of task
        % spikeStruct.ts is already sorted
        numSpikesTotal = numel(spikeStruct.ts);
        pctile25 = round(0.25 * numSpikesTotal);
        pctile75 = round(0.75 * numSpikesTotal);
        numSpikesMain = pctile75 - pctile25;
        totalTimeMain = spikeStruct.ts(pctile75) - spikeStruct.ts(pctile25);
        meanFiringRates(unitCount) = numSpikesMain / totalTimeMain;
        
        if meanFiringRates(unitCount) < minFiringRate
            isCell(unitCount) = false;
            isBroadSpiking(unitCount) = false;
            isNarrowSpiking(unitCount) = false;
        end
    end
end

%% plot
fprintf('-------------------------------\n');
subdivisions = {'PLd', 'PM', 'PI', 'PUL', 'all'};
for i = 1:numel(subdivisions)
    subdivision = subdivisions{i};
    if strcmp(subdivision, 'all')
        isInSubdivision = true(size(troughToPeakTimes));
    elseif strcmp(subdivision, 'PUL')
        isInSubdivision = strcmp(localization, 'PLd') | strcmp(localization, 'PLv') | ...
                strcmp(localization, 'PM') | strcmp(localization, 'PI');
    else
        isInSubdivision = strcmp(localization, subdivision);
    end
    
    numBS = sum(isBroadSpiking & isInSubdivision);
    numNS = sum(isNarrowSpiking & isInSubdivision);
    numCells = numBS + numNS;
    assert(numCells == sum(isCell & isInSubdivision));
    
    % deal with rounding issues. data is at intervals of 0.025 ms
    troughToPeakSamples = round(troughToPeakTimes * 40);
    
    figure_tr_inch(8, 6); clf;
    subaxis(1, 1, 1, 'ML', 0.08, 'MR', 0.03);
    set(gcf, 'Color', 'white');
    set(gcf, 'renderer', 'painters');
    hold on;
    binEdges = (0:0.025:0.7) + 0.025/2;
    binEdgesSamples = binEdges * 40;
    histogram(troughToPeakSamples(isBroadSpiking & isInSubdivision), binEdgesSamples);
    histogram(troughToPeakSamples(isNarrowSpiking & isInSubdivision), binEdgesSamples);
    xTickSamples = (0.1:0.1:0.7) * 40;
    set(gca, 'XTick', xTickSamples);
    set(gca, 'XTickLabel', xTickSamples / 40);
    ylabel('Number of Cells');
    xlabel('Trough to Peak Time (ms)');
    box off;
%     grid on;
    set(gca, 'FontSize', 12);
    title(sprintf('%s (N=%d)', subdivision, numCells), 'FontSize', 16);

    plotFileName = sprintf('%s/allSessions-%s-troughToPeakTimeHist.png', processedDataRootDir, subdivision);
    export_fig(plotFileName, '-nocrop');
    % close;

    fprintf('%s: There are %d/%d = %d%% BS cells and %d/%d = %d%% NS cells.\n', subdivision, ...
            numBS, numCells, round(numBS / numCells * 100), ...
            numNS, numCells, round(numNS / numCells * 100));
    fprintf('%s: Max trough to peak time is %0.3f ms.\n', subdivision, max(troughToPeakTimes(isCell & isInSubdivision)));
end

%% plot
fprintf('-------------------------------\n');
subdivisions = {'PLd', 'PM', 'PI', 'PUL', 'all'};
for i = 1:numel(subdivisions)
    subdivision = subdivisions{i};
    if strcmp(subdivision, 'all')
        isInSubdivision = true(size(meanFiringRates));
    elseif strcmp(subdivision, 'PUL')
        isInSubdivision = strcmp(localization, 'PLd') | strcmp(localization, 'PLv') | ...
                strcmp(localization, 'PM') | strcmp(localization, 'PI');
    else
        isInSubdivision = strcmp(localization, subdivision);
    end
    
    numBS = sum(isBroadSpiking & isInSubdivision);
    numNS = sum(isNarrowSpiking & isInSubdivision);
    numCells = numBS + numNS;
    assert(numCells == sum(isCell & isInSubdivision));
    
    figure_tr_inch(8, 6); clf;
    subaxis(1, 1, 1, 'ML', 0.08, 'MR', 0.03);
    set(gcf, 'Color', 'white');
    set(gcf, 'renderer', 'painters');
    hold on;
    binEdges = 0:2:60;
    histogram(meanFiringRates(isBroadSpiking & isInSubdivision), binEdges);
    histogram(meanFiringRates(isNarrowSpiking & isInSubdivision), binEdges);
    ylabel('Number of Cells');
    xlabel('Mean Firing Rate (Hz)');
    box off;
%     grid on;
    set(gca, 'FontSize', 12);
    title(sprintf('%s (N=%d)', subdivision, numCells), 'FontSize', 16);

    plotFileName = sprintf('%s/allSessions-%s-meanFiringRateHist.png', processedDataRootDir, subdivision);
    export_fig(plotFileName, '-nocrop');
    % close;

    fprintf('%s: There are %d/%d = %d%% BS cells and %d/%d = %d%% NS cells.\n', subdivision, ...
            numBS, numCells, round(numBS / numCells * 100), ...
            numNS, numCells, round(numNS / numCells * 100));
    fprintf('%s: Max firing rate is %0.1f Hz.\n', subdivision, max(meanFiringRates(isCell & isInSubdivision)));
    fprintf('%s: Median firing rate for BS cells (N=%d) is %0.1f Hz.\n', subdivision, ...
            numBS, median(meanFiringRates(isBroadSpiking & isInSubdivision)));
    fprintf('%s: Median firing rate for NS cells (N=%d) is %0.1f Hz.\n', subdivision, ...
            numNS, median(meanFiringRates(isNarrowSpiking & isInSubdivision)));
%     ranksum(meanFiringRates(isBroadSpiking & isInSubdivision), meanFiringRates(isNarrowSpiking & isInSubdivision))
end

%%
fprintf('-------------------------------\n');
subdivisions = {'PUL'};
for i = 1:numel(subdivisions)
    subdivision = subdivisions{i};
    if strcmp(subdivision, 'all')
        isInSubdivision = true(size(meanFiringRates));
    elseif strcmp(subdivision, 'PUL')
        isInSubdivision = strcmp(localization, 'PLd') | strcmp(localization, 'PLv') | ...
                strcmp(localization, 'PM') | strcmp(localization, 'PI');
    else
        isInSubdivision = strcmp(localization, subdivision);
    end
    
    numBS = sum(isBroadSpiking & isInSubdivision);
    numNS = sum(isNarrowSpiking & isInSubdivision);
    numCells = numBS + numNS;
    assert(numCells == sum(isCell & isInSubdivision));
    
    figure_tr_inch(5, 4); clf;
    subaxis(1, 1, 1, 'ML', 0.13, 'MR', 0.03);
    set(gcf, 'Color', 'white');
    set(gcf, 'renderer', 'painters');
    hold on;
    cols = lines(2);
    nWfTime = size(meanWfs, 2);
    spikeFs = 40000;
    waveformT = (1:nWfTime)/(spikeFs/1000);
    meanWfsBS = meanWfs(isBroadSpiking & isInSubdivision,:);
    meanWfsNS = meanWfs(isNarrowSpiking & isInSubdivision,:);
    for j = 1:size(meanWfsBS, 1)
        [~,troughInd] = min(meanWfsBS(j,:));
        plot(waveformT - troughInd/(spikeFs/1000), meanWfsBS(j,:), 'Color', cols(1,:));
    end
    for j = 1:size(meanWfsNS, 1)
        [~,troughInd] = min(meanWfsNS(j,:));
        plot(waveformT - troughInd/(spikeFs/1000), meanWfsNS(j,:), 'Color', cols(2,:));
    end
    ylabel('Voltage (mV)');
    xlabel('Time from Trough (ms)');
    xlim([-0.3 0.9]);
    ylim([-0.12 0.12]);
    box off;
%     grid on;
    set(gca, 'FontSize', 12);
    title(sprintf('%s (N=%d)', subdivision, numCells), 'FontSize', 16);

    plotFileName = sprintf('%s/allSessions-%s-meanWaveforms.png', processedDataRootDir, subdivision);
    export_fig(plotFileName, '-nocrop');
    % close;

end