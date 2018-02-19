function muaVEPMapping(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad)
% MUA VEP Mapping, one channel

%% setup and load data
v = 10;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('VEP Mapping Analysis - MUA\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('MUA Channel to Load: %d\n', muaChannelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

doPlot = 1;

%% input check
assert(numel(muaChannelsToLoad) == 1);

%% load recording information
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, muaChannelsToLoad, 'VEPM', 1, 0);
sessionName = R.sessionName;

%% find nonsparse blocks
isFiringNonsparseByBlock = findNonsparseBlocks(D, D.allMUAStructs, R.vepmIndices);

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
else
    origFlashEvents = D.events{3};
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

numRandomizations = 500;

%% iterate through each unit
nUnits = numel(D.allMUAStructs);
assert(nUnits == 1);
j = 1; % just one unit at a time in this script

%% compute evoked spiking and make SPDF plots for visual and motor events
fprintf('-------------------------------------------------------------\n');
fprintf('Processing %d MUAs...\n', nUnits);

muaStruct = D.allMUAStructs{j};
unitName = muaStruct.name;
spikeTimes = muaStruct.ts;
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
    blockStartTime = D.blockStartTimes(R.vepmIndices(sparseBlocks(k)));
    blockStopTime = D.blockStopTimes(R.vepmIndices(sparseBlocks(k)));
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
vepmIndices = R.vepmIndices;
muaStruct.vepmPsthParams = var2struct(vepmIndices, psthWindow, ...
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
    
%% plot
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
    plotSpikeWaveform(D, j, 1);

    %% info
    writeUnitInfo(muaStruct, axBig, -0.03, infoTextTop);

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
muaStruct.isVisuallyDriven = -1;
muaStruct.peakLatencyInfo = NaN;
muaStruct.troughLatencyInfo = NaN;
muaStruct.latencyInfo = NaN;

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
        muaStruct.isVisuallyDriven = 1;
    else
        visRespText = 'Not Visually Driven';
        visRespTextCol = [0.9 0.2 0.05];
        muaStruct.isVisuallyDriven = 0;
    end
    
    %% compute latency
    minPeakForLatency = minAbsNormMaxRespVisuallyDriven;
    maxTroughForLatency = -minAbsNormMaxRespVisuallyDriven;
    latencyCol = zeros(3, 1);
    
    if muaStruct.isVisuallyDriven
        % check peak
        peakLatencyInfo = computeLatencyPeakMethod(normPsthResponse, psthT, ...
                psthWindow, latencyIndices, 0, minPeakForLatency, 0);
        
        % check trough
        troughLatencyInfo = computeLatencyPeakMethod(normPsthResponse, psthT, ...
                psthWindow, latencyIndices, 0, maxTroughForLatency, 1);
        
        % get the earlier one
        latencyInfo = getEarlierPeakTroughLatencyInfo(peakLatencyInfo, troughLatencyInfo);
        
        muaStruct.peakLatencyInfo = peakLatencyInfo;
        muaStruct.troughLatencyInfo = troughLatencyInfo;
        muaStruct.latencyInfo = latencyInfo;
        
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
    
    fprintf('\t\t%s - %s, %s, Latency = %s\n', muaStruct.name, ...
            muaStruct.physClass, visRespText, latencyText);

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
    fprintf('\t\t%s - No response.\n', muaStruct.name);
    
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

%% save plot
if doPlot
    plotFileName = sprintf('%s/%s-%s-vepm-v%d.png', processedDataDir, unitName, blockName, v);
    fprintf('Saving to %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
    close;
end

%% save data
saveFileName = sprintf('%s/%s-%s-vepm-v%d.mat', processedDataDir, unitName, blockName, v);
fprintf('Saving to %s...\n', saveFileName);
save(saveFileName, 'muaStruct');

