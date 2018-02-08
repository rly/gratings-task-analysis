function spikeVEPMapping(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, spikeChannelsToLoad)

%% setup and load data

v = 10;

fprintf('\n-------------------------------------------------------\n');
fprintf('VEP Mapping Analysis - Spikes\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('Spike Channels to Load: %d\n', spikeChannelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

doPlot = 1;

%% load recording information
[R, D, sessionName, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad);

%% find nonsparse blocks
isFiringNonsparseByBlock = findNonsparseBlocks(D, R.vepmIndices);

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

%% SPDF params

psthWindow = [0.25 0.3]; % secs before, secs after

analysisWindowOffset = [0.025 0.225]; % secs offset from event onset
analysisWindow = analysisWindowOffset + psthWindow(1);
analysisWindowDur = diff(analysisWindowOffset);

baselineWindowOffset = [-0.2 0];
baselineWindow = baselineWindowOffset + psthWindow(1);
baselineWindowDur = diff(baselineWindowOffset);

latencyWindowOffset = [0.025 0.125]; % secs offset from event onset
latencyWindow = latencyWindowOffset + psthWindow(1);

kernelSigma = 0.01;
% nTime = fix(5*sum(psthWindow)/kernelSigma) + 1; % number of time steps. 
nTime = fix(10*sum(psthWindow)/kernelSigma) + 1; % number of time steps. 
psthT = linspace(0, sum(psthWindow), nTime); % starts at 0
t = psthT - psthWindow(1);

analysisIndices = getTimeLogicalWithTolerance(psthT, analysisWindow);
baselineIndices = getTimeLogicalWithTolerance(psthT, baselineWindow);
latencyIndices = getTimeLogicalWithTolerance(psthT, latencyWindow);

minAbsNormMaxRespVisuallyDriven = 3; % normed value is z-scored (SD units)
minProportionFlashesWithSpikeVisuallyDriven = 0.1;

%% iterate through each unit
nUnits = numel(D.allSpikeStructs);
fprintf('Processing %d units...\n', nUnits);

for j = 1:nUnits
    spikeStruct = D.allSpikeStructs{j};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
            nUnits, round(j/nUnits*100));

%% compute psth for each flash condition
% this does not subtract out per-condition baseline activity

% remove spike and event times during sparse blocks
spikeTimesClean = spikeTimes;
origFlashEventsClean = origFlashEvents;
sparseBlocks = find(~isFiringNonsparseByBlock(:,j));
rangeSparseRemovedFlashEvents = nan(numel(sparseBlocks), 2); % start, end
for k = 1:numel(sparseBlocks)
    blockStartTime = D.blockStartTimes(vepmIndices(sparseBlocks(k)));
    blockStopTime = D.blockStopTimes(vepmIndices(sparseBlocks(k)));
    spikeTimesClean(spikeTimesClean >= blockStartTime & spikeTimesClean <= blockStopTime) = [];
    origFlashEventsClean(origFlashEventsClean >= blockStartTime & origFlashEventsClean <= blockStopTime) = [];
    rangeSparseRemovedFlashEvents(k,:) = [find(origFlashEvents >= blockStartTime, 1, 'first') find(origFlashEvents <= blockStopTime, 1, 'last')];
end

allAlignedSpikeTimes = createdatamatpt(spikeTimes, origFlashEvents, psthWindow);
allAlignedSpikeTimesClean = createdatamatpt(spikeTimesClean, origFlashEventsClean, psthWindow);

% count spikes in analysis and baseline windows
nSpikesAnalysis = zeros(numel(allAlignedSpikeTimesClean), 1);
nSpikesBaseline = zeros(numel(allAlignedSpikeTimesClean), 1);
for n = 1:numel(allAlignedSpikeTimesClean)
    nSpikesAnalysis(n) = sum(analysisWindow(1) <= allAlignedSpikeTimesClean(n).times & ...
        allAlignedSpikeTimesClean(n).times <= analysisWindow(2));
    nSpikesBaseline(n) = sum(baselineWindow(1) <= allAlignedSpikeTimesClean(n).times & ...
        allAlignedSpikeTimesClean(n).times <= baselineWindow(2));
end
proportionFlashesWithSpike = sum(nSpikesAnalysis > 0) / numel(nSpikesAnalysis);

% if the cell is quiet for 20+ consecutive flashes in both analysis and
% baseline windows, remove those flashes
% consecSilentTrialsToRemove = 20;
% totalSpikesAnalyzedPerTrial = nSpikesAnalysis + nSpikesBaseline;
% consecSilentTrials = 0;
% isRemoveTrials = false(numel(totalSpikesAnalyzedPerTrial), 1);
% for n = 1:numel(totalSpikesAnalyzedPerTrial)
%     if totalSpikesAnalyzedPerTrial(n) == 0
%         consecSilentTrials = consecSilentTrials + 1;
%     else
%         % when a non-silent trial is encountered, then remove the preceding
%         % silent trials if there are enough of them
%         if consecSilentTrials >= consecSilentTrialsToRemove
%             isRemoveTrials(n-consecSilentTrials:n-1) = 1;
%         end
%         consecSilentTrials = 0;
%     end
% end
% if consecSilentTrials >= consecSilentTrialsToRemove
%     isRemoveTrials(n-consecSilentTrials+1:n) = 1;
% end
% nSpikesAnalysis(isRemoveTrials) = [];
% nSpikesBaseline(isRemoveTrials) = [];
% allAlignedSpikeTimesCleaner = allAlignedSpikeTimesClean;
% allAlignedSpikeTimesCleaner(isRemoveTrials) = [];

% compute metrics using spike rate by COUNTing spikes and dividing by time
baselineMetricByCount = mean(nSpikesBaseline) / baselineWindowDur; % average num spikes per trial
baselineMetricSDOverTrialsByCount = std(nSpikesBaseline) / baselineWindowDur;
baselineMetricSEMByCount = baselineMetricSDOverTrialsByCount / sqrt(size(nSpikesAnalysis, 1));
responseMetricByCount = mean(nSpikesAnalysis) / analysisWindowDur; % average num spikes per trial
responseMetricSDOverTrialsByCount = std(nSpikesAnalysis) / analysisWindowDur;

% stats: what is the probability of observing this number of spikes in the
% analysis window if we assume the number of spikes in the analysis window
% and in the baseline window come from the same distribution?
[~,responseMetricPValueByCount,~,responseMetricStatsByCount] = ttest(nSpikesAnalysis / analysisWindowDur, nSpikesBaseline / baselineWindowDur);

% compute how many SDs away is the response
if baselineMetricSDOverTrialsByCount > 1e-8
    normResponseMetricByCount = (responseMetricByCount - baselineMetricByCount) / baselineMetricSDOverTrialsByCount;
else
    normResponseMetricByCount = 0;
end

% compute psth response
[psthResponse,~,psthBootstrapErr] = fixedPsth(allAlignedSpikeTimesClean, kernelSigma, 2, psthT);

% compute metrics using PSTH estimated spike rate
if ~all(isnan(psthResponse))
    baselineResponseByPsth = mean(psthResponse(baselineIndices));
    % SEMByPsth & SDOverTimeByPsth are correlated but not equal
    % SDOverTime makes more sense when seeing whether a time point in the
    % average deviates from the norm
    baselineResponseSEMByPsth = mean(psthBootstrapErr(baselineIndices));
    baselineResponseSDOverTimeByPsth = std(psthResponse(baselineIndices));
    
    meanResponseByPsth = mean(psthResponse(analysisIndices));
    maxResponseByPsth = max(psthResponse(analysisIndices));
    minResponseByPsth = min(psthResponse(analysisIndices));
    
    % compute how many SDs away is the mean/max/min response
    if baselineResponseSDOverTimeByPsth > 1e-8
        normPsthResponse = (psthResponse - baselineResponseByPsth) / baselineResponseSDOverTimeByPsth;
        normMeanResponseByPsth = (meanResponseByPsth - baselineResponseByPsth) / baselineResponseSDOverTimeByPsth;
        normMaxResponseByPsth = (maxResponseByPsth - baselineResponseByPsth) / baselineResponseSDOverTimeByPsth;
        normMinResponseByPsth = (minResponseByPsth - baselineResponseByPsth) / baselineResponseSDOverTimeByPsth;
        assert(abs(normMeanResponseByPsth - mean(normPsthResponse(analysisIndices))) < 1e-5); % double check
    else
        normPsthResponse = 0;
        normMeanResponseByPsth = 0;
        normMaxResponseByPsth = 0;
        normMinResponseByPsth = 0;
    end
    
    % compute non-parametric stats on psth response
    % bootstrap on mean of baseline PSTH with 500 shuffles
    % compare actual mean psth response to distribution of bootstrapped
    % baselines
    numRandomizations = 200;
    bootstrappedMeanBaselines = zeros(numRandomizations, 1);
    numTrialsClean = numel(allAlignedSpikeTimesClean);
    for m = 1:numRandomizations
        bootRandInd = randi(numTrialsClean, 1, numTrialsClean);
        bootstrapPsthResponse = fixedPsth(allAlignedSpikeTimesClean(bootRandInd), kernelSigma, 0, psthT);
        if ~isempty(bootstrapPsthResponse)
            bootstrappedMeanBaselines(m) = mean(bootstrapPsthResponse(baselineIndices));
        end % else, then there were no spikes and leave as zeros
    end
    if meanResponseByPsth > baselineResponseByPsth
        responseMetricPValueByBootstrapPsth = sum(meanResponseByPsth < bootstrappedMeanBaselines) / numRandomizations * 2;
    else
        responseMetricPValueByBootstrapPsth = sum(meanResponseByPsth > bootstrappedMeanBaselines) / numRandomizations * 2;
    end
else
    % nan out all these vars...
    nantemp = num2cell(nan(13, 1));
    [baselineResponseByPsth, ...
        baselineResponseSEMByPsth, baselineResponseSDOverTimeByPsth, ...
        meanResponseByPsth, normPsthResponse, normMeanResponseByPsth, ...
        maxResponseByPsth, normMaxResponseByPsth, ...
        minResponseByPsth, normMinResponseByPsth, ...
        numRandomizations, bootstrappedMeanBaselines, ...
        responseMetricPValueByBootstrapPsth] = deal(nantemp{:});
end

%% save ALL the values
% there's gotta be a better way
D.allSpikeStructs{j}.vepmPsthParams = var2struct(vepmIndices, psthWindow, ...
        analysisWindowOffset, baselineWindowOffset, ...
        latencyWindowOffset, kernelSigma, psthT, t, spikeTimesClean, ...
        origFlashEventsClean, allAlignedSpikeTimes, allAlignedSpikeTimesClean, ...
        baselineMetricByCount, baselineMetricSDOverTrialsByCount, ...
        baselineMetricSEMByCount, responseMetricByCount, ...
        responseMetricSDOverTrialsByCount, responseMetricPValueByCount, ...
        responseMetricStatsByCount, normResponseMetricByCount, ...
        psthResponse, psthBootstrapErr, baselineResponseByPsth, ...
        baselineResponseSEMByPsth, baselineResponseSDOverTimeByPsth, ...
        meanResponseByPsth, normPsthResponse, normMeanResponseByPsth, ...
        maxResponseByPsth, normMaxResponseByPsth, ...
        minResponseByPsth, normMinResponseByPsth, ...
        numRandomizations, bootstrappedMeanBaselines, ...
        responseMetricPValueByBootstrapPsth, minAbsNormMaxRespVisuallyDriven);
    
%%
if doPlot
    f = figure_tr_inch(13, 7.5); clf;
    set(gcf, 'Color', 'white');
    set(gcf, 'renderer', 'painters');

    %% make main title
    axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
    set(get(axBig, 'Title'), 'Visible', 'on')

    modTitle = sprintf('Flash VES: %s (%d flashes)', unitName, nFlashes);
    titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
    title(modTitle, 'FontSize', 14, titleParams{:});

    %% location params
    waveformW = 0.1;
    rasterByTimeW = 0.75;
    spdfW = 0.75;
    waveformH = 0.15;
    rasterByTimeH = 0.41;
    spdfH = 0.41;

    waveformLeft1 = 0.05;
    rasterByTimeLeft1 = waveformLeft1 + waveformW + 0.08;
    spdfLeft = rasterByTimeLeft1;
    btm = 0.07;
    infoTextTop = btm + 0.1;
    infoText2Top = infoTextTop + 0.53;
    rasterByTimeBtm = btm + spdfH + 0.04;
    waveformBtm = 0.76;

    %% plot spike waveform
    axes('Position', [waveformLeft1 waveformBtm waveformW waveformH]); 
    plotSpikeWaveform(D, j);

    %% info
    writeUnitInfo(spikeStruct, axBig, -0.03, infoTextTop);

    %% plot raster
    axes('Position', [rasterByTimeLeft1 rasterByTimeBtm rasterByTimeW rasterByTimeH]); 
    hold on;

    xBounds = [-0.2 0.27];
    window = psthWindow;
    data = allAlignedSpikeTimes;
    rasterY = 0;
    lineParams = {'Color', [0 0 0], 'LineWidth', 1};
    lineHeight = 5;
    for k = 1:numel(data)
        rasterY = rasterY + 1;
        if ~isempty(data(k).times)
            adjTimes = data(k).times - window(1);
            for l = 1:numel(data(k).times)
                % plot mini lines at each spike time
                plot([adjTimes(l) adjTimes(l)], rasterY + [-lineHeight lineHeight]/2, lineParams{:});
            end
        end
        % plot red asterisks at removed sets of trials due to no firing
%         if isRemoveTrials(k)
%             plot(xBounds(1), rasterY, 'r*');
%         end
    end

    % plot event line marker
    yBounds = [0 rasterY+1];
    plot([0 0], yBounds, '-', 'Color', 0.5*ones(3, 1));
    
    % shade in red removed sparse blocks
    for k = 1:numel(sparseBlocks)
        fillH = jbfill(xBounds, rangeSparseRemovedFlashEvents(k,1) + [-0.5 -0.5], ...
                rangeSparseRemovedFlashEvents(k,2) + [0.5 0.5], ...
                [1 0 0], ones(3, 1), 0.1);
        uistack(fillH, 'bottom');
    end
    
    % shade in baseline and analysis windows
    fillH = jbfill(baselineWindowOffset, [0 0], [yBounds(2) yBounds(2)], [0.7 0.7 1], ones(3, 1), 0.1);
    uistack(fillH, 'bottom');
    fillH = jbfill(analysisWindowOffset, [0 0], [yBounds(2) yBounds(2)], [0.3 1 0.3], ones(3, 1), 0.1);
    uistack(fillH, 'bottom');
    
    xlim(xBounds);
    ylim(yBounds);
    set(gca, 'YDir', 'reverse');
    set(gca, 'XTickLabel', []);
    ylabel('Flash number');
    set(gca, 'Layer', 'top');

    %% spdf
    spdfAx = axes('Position', [spdfLeft btm spdfW spdfH]); 
    if ~isempty(psthResponse)
        % plot baseline and x SDs around baseline
        hold on;
        plot(xBounds, baselineResponseByPsth * ones(size(xBounds)), 'm');
        baselineShadingUB = baselineResponseByPsth + minAbsNormMaxRespVisuallyDriven * baselineResponseSDOverTimeByPsth;
        baselineShadingLB = baselineResponseByPsth - minAbsNormMaxRespVisuallyDriven * baselineResponseSDOverTimeByPsth;
        fillH = jbfill(xBounds, baselineShadingUB * ones(size(xBounds)), ...
                baselineShadingLB * ones(size(xBounds)), [1 0 1], ones(3, 1), 0.08);
        uistack(fillH, 'bottom');

        % plot SPDF and bootstrap error
        hold on;
        plot(t, psthResponse, 'Color', lines(1), 'LineWidth', 5);
        shadingUB = psthResponse + psthBootstrapErr;
        shadingLB = psthResponse - psthBootstrapErr;
        fillH = jbfill(t, shadingUB, shadingLB, lines(1), ones(3, 1), 0.2);
        uistack(fillH, 'bottom');

        % plot event line marker
        hold on;
        origYLim = ylim();
        yBounds = [min(0, origYLim(1)) max(1, origYLim(2))];
        plot([0 0], yBounds, '-', 'Color', 0.5*ones(3, 1)); 

        % shade in baseline and analysis windows
        fillH = jbfill(baselineWindowOffset, [yBounds(1) yBounds(1)], ...
                [yBounds(2) yBounds(2)], [0.7 0.7 1], ones(3, 1), 0.1);
        uistack(fillH, 'bottom');
        fillH = jbfill(analysisWindowOffset, [yBounds(1) yBounds(1)], ...
                [yBounds(2) yBounds(2)], [0.3 1 0.3], ones(3, 1), 0.1);
        uistack(fillH, 'bottom');
        hold on;

        xlim(xBounds);
        ylim(yBounds);
        xlabel('Time from flash onset (s)');
        ylabel('Estimated spike rate (Hz)');
    end
end

%% add response metrics text
D.allSpikeStructs{j}.isVisuallyDriven = -1;
D.allSpikeStructs{j}.peakLatencyInfo = NaN;
D.allSpikeStructs{j}.troughLatencyInfo = NaN;
D.allSpikeStructs{j}.latencyInfo = NaN;

if ~isempty(psthResponse)
    %% classify cell
    % TODO use a proper statistical test to see if there is a visual
    % response
    % note that units may be visually responsive and we can't tell because 
    % there is not enough firing during this task
    
    if proportionFlashesWithSpike >= minProportionFlashesWithSpikeVisuallyDriven && ...
            (abs(normMaxResponseByPsth) >= minAbsNormMaxRespVisuallyDriven || ...
            abs(normMinResponseByPsth) >= minAbsNormMaxRespVisuallyDriven)
        visRespText = 'Visually Driven';
        visRespTextCol = [0.05 0.5 0.05];
        D.allSpikeStructs{j}.isVisuallyDriven = 1;
    else
        visRespText = 'Not Visually Driven';
        visRespTextCol = [0.9 0.2 0.05];
        D.allSpikeStructs{j}.isVisuallyDriven = 0;
    end
    
    %% compute latency
    minPeakForLatency = minAbsNormMaxRespVisuallyDriven;
    maxTroughForLatency = -minAbsNormMaxRespVisuallyDriven;
    latencyCol = zeros(3, 1);
    
    if D.allSpikeStructs{j}.isVisuallyDriven
        % check peak
        peakLatencyInfo = computeLatencyPeakMethod(normPsthResponse, psthT, ...
                psthWindow, latencyIndices, 0, minPeakForLatency, 0);
        
        % check trough
        troughLatencyInfo = computeLatencyPeakMethod(normPsthResponse, psthT, ...
                psthWindow, latencyIndices, 0, maxTroughForLatency, 1);
        
        % get the earlier one
        latencyInfo = getEarlierPeakTroughLatencyInfo(peakLatencyInfo, troughLatencyInfo);
        
        D.allSpikeStructs{j}.peakLatencyInfo = peakLatencyInfo;
        D.allSpikeStructs{j}.troughLatencyInfo = troughLatencyInfo;
        D.allSpikeStructs{j}.latencyInfo = latencyInfo;
        
        if doPlot && ~isnan(latencyInfo.latency)
            latencyCol = [0.6 0.05 0.6];
            plot(latencyInfo.latency, psthResponse(latencyInfo.latencyTInd), 'o', ...
                    'Color', latencyCol, 'MarkerSize', 8, 'LineWidth', 2);
            % NOTE that the latency may be a little later than the exact 
            % half peak depending on the time resolution
        end
        latencyText = sprintf('%d ms', round(latencyInfo.latency * 1000));
    else
        latencyText = 'None';
    end
    
    fprintf('\t\t%s - %s, %s, Latency = %s\n', D.allSpikeStructs{j}.name, ...
            D.allSpikeStructs{j}.physClass, visRespText, latencyText);

    %% text
    if doPlot
        if responseMetricPValueByBootstrapPsth < 0.05
            responseMetricPValueByBootstrapPsthCol = [0.05 0.5 0.05];
        else
            responseMetricPValueByBootstrapPsthCol = [0.9 0.2 0.05];
        end
        
        textParams = {'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top'};
        text(axBig, -0.03, infoText2Top, {...
                sprintf('Blocks: %s', blockName), ...
                '', ...
                ...sprintf('Baseline (count): %0.1f Hz', baselineMetricByCount), ...
                ...sprintf('Response (count): %0.1f Hz (p = %0.3f)', responseMetricByCount, responseMetricPValueByCount), ...
                ...'', ...
                sprintf('Flashes w/ Spike Response = %d%%', round(proportionFlashesWithSpike * 100)), ...
                sprintf('Mean Response (SPDF) = %0.1f Hz', meanResponseByPsth), ...
                sprintf('Mean Baseline (SPDF) = %0.1f Hz', baselineResponseByPsth), ...
                sprintf('Response > Baseline Boot. {\\color[rgb]{%f,%f,%f}P = %0.3f}', responseMetricPValueByBootstrapPsthCol, responseMetricPValueByBootstrapPsth), ...
                '', ...
                sprintf('SD Baseline Over Time = %0.1f Hz', baselineResponseSDOverTimeByPsth), ...
                sprintf('Max Response = %0.1f Hz (%0.1f SDs)', maxResponseByPsth, normMaxResponseByPsth), ...
                sprintf('Min Response = %0.1f Hz (%0.1f SDs)', minResponseByPsth, normMinResponseByPsth), ...
                '', ...
                sprintf('"Visually Driven": \\geq %d%% Flashes w/ Spike', round(minProportionFlashesWithSpikeVisuallyDriven * 100)), ...
                sprintf('  and Max or Min Response \\geq %0.1f SDs', minAbsNormMaxRespVisuallyDriven), ...
                sprintf('Functional Class: {\\color[rgb]{%f,%f,%f}%s}', visRespTextCol, visRespText), ...
                '', ...
                sprintf('Latency Window: [%0.3f, %0.3f] s', latencyWindowOffset), ...
                sprintf('Latency Peak/Trough \\geq %0.1f SDs', minPeakForLatency), ...
                sprintf('Latency (Time to Half-Peak) = {\\color[rgb]{%f,%f,%f}%s}', latencyCol, latencyText)}, ...
                textParams{:});
    end
else
    fprintf('\t\t%s - No response.\n', D.allSpikeStructs{j}.name);
    
    %% text
    if doPlot
        textParams = {'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top'};
        text(axBig, -0.03, infoText2Top, {...
                sprintf('Blocks: %s', blockName), ...
                '', ...
                sprintf('Baseline Window: [%0.2f, %0.2f] ms', baselineWindowOffset), ...
                sprintf('Response Window: [%0.2f, %0.2f] ms', analysisWindowOffset)}, ...
                textParams{:});
    end
end

%% save
if doPlot
    plotFileName = sprintf('%s/%s-%s.png', processedDataDir, unitName, blockName);
    export_fig(plotFileName, '-nocrop');
    close;
end

clear spikeTimesClean ...
        origFlashEventsClean allAlignedSpikeTimes allAlignedSpikeTimesClean ...
        baselineMetricByCount baselineMetricSDOverTrialsByCount ...
        baselineMetricSEMByCount responseMetricByCount ...
        responseMetricSDOverTrialsByCount responseMetricPValueByCount ...
        responseMetricStatsByCount normResponseMetricByCount ...
        psthResponse psthBootstrapErr baselineResponseByPsth ...
        baselineResponseSEMByPsth baselineResponseSDOverTimeByPsth ...
        meanResponseByPsth normPsthResponse normMeanResponseByPsth ...
        numRandomizations bootstrappedMeanBaselines ...
        responseMetricPValueByBootstrapPsth;

end % for each unit

%% write updated D struct
% remove wf and ts to save space
for j = 1:nUnits
    D.allSpikeStructs{j}.wf = [];
    D.allSpikeStructs{j}.ts = [];
end
saveFileName = sprintf('%s/D-afterSpikeVEPM-%s.mat', processedDataDir, blockName);
fprintf('Writing post-processed data file: %s ...', saveFileName);
tic;
save(saveFileName, 'D');
fprintf(' done (%0.2f s).\n', toc);

%% plot cdf of latencies
nUnits = numel(D.allSpikeStructs);
fprintf('Processing %d units...\n', nUnits);

latencies = nan(nUnits, 3); % latency, channel number, spike struct index
for j = 1:nUnits
    spikeStruct = D.allSpikeStructs{j};
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
meanLatenciesByChannel = nan(D.nSpikeCh, 1);
for i = 1:numel(meanLatenciesByChannel)
    meanLatenciesByChannel(i) = nanmean(latencies(latencies(:,2) == spikeChannelsToLoad(i)), 1);
end

channelIndices = spikeChannelsToLoad;
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
    spikeStruct = D.allSpikeStructs{j};
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
    spikeStruct = D.allSpikeStructs{j};
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
    spikeStruct = D.allSpikeStructs{j};
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
    spikeStruct = D.allSpikeStructs{j};
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


%%
bsLatencies = nan(nUnits, 1);
nsLatencies = nan(nUnits, 1);

for j = 1:nUnits
    spikeStruct = D.allSpikeStructs{j};
    
    % include all units even if there was insignificant response
    % but exclude if no firing or super sparse firing
    if any(isnan(spikeStruct.vepmPsthParams.normPsthResponse)) || ...
            spikeStruct.vepmPsthParams.meanResponseByPsth < minResponseByPsth
        continue;
    end
    
    if ~isstruct(spikeStruct.latencyInfo)
        continue;
    end
    
    % exclude if axon
    if ~(strcmp(spikeStruct.physClass, 'Broad-Spiking') || strcmp(spikeStruct.physClass, 'Narrow-Spiking'))
        continue;
    end
    
    if strcmp(spikeStruct.physClass, 'Broad-Spiking')
        bsLatencies(j) = spikeStruct.latencyInfo.latency * 1000;
    else
        nsLatencies(j) = spikeStruct.latencyInfo.latency * 1000;
    end
end
bsLatencies(isnan(bsLatencies)) = [];
nsLatencies(isnan(nsLatencies)) = [];

fprintf('\n');
fprintf('BS Cells (N = %d): Mean = %d ms, Median = %d ms\n', numel(bsLatencies), round(mean(bsLatencies)), round(median(bsLatencies)));
fprintf('NS Cells (N = %d): Mean = %d ms, Median = %d ms\n', numel(nsLatencies), round(mean(nsLatencies)), round(median(nsLatencies)));
fprintf('Excluded Units: %d\n', nUnits - numel(bsLatencies) - numel(nsLatencies));

%% write update D struct for Leenoy
% tic;
% isLoadSpikes = 1;
% isLoadLfp = 0;
% isLoadSpkc = 0;
% D2 = loadPL2(pl2FileName, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc); 
% for j = 1:nUnits
%     D2.allSpikeStructs{j}.wf = [];
%     D2.allSpikeStructs{j}.isVisuallyDriven = D.allSpikeStructs{j}.isVisuallyDriven;
%     D2.allSpikeStructs{j}.peakLatencyInfo = D.allSpikeStructs{j}.peakLatencyInfo;
%     D2.allSpikeStructs{j}.troughLatencyInfo = D.allSpikeStructs{j}.troughLatencyInfo;
%     D2.allSpikeStructs{j}.latencyInfo = D.allSpikeStructs{j}.latencyInfo;
%     D2.allSpikeStructs{j}.vepmPsthParams = D.allSpikeStructs{j}.vepmPsthParams;
% end
% D = D2;
% saveFileName = sprintf('%s/D-all-withVisuallyDrivenTag.mat', processedDataDir);
% fprintf('Writing post-processed data file: %s ...', saveFileName);
% tic;
% save(saveFileName, 'D');
% fprintf(' done (%0.2f s).\n', toc);