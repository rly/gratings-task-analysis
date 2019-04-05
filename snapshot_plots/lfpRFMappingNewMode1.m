function lfpRFMappingNewMode1(processedDataRootDir, dataDirRoot, muaDataDirRoot, ...
        recordingInfoFileName, sessionInd, channelsToLoad, rfMappingNewInfoFileName)
% LFP RF Mapping, all channels on a probe
% can't really do one channel at a time because of Common Average
% Referencing

% TODO: incorporate shuffle test. what is probability of a peak in smoothed
% map of that amplitude if there is no difference between conditions
% computationally intensive: shuffle trials
% easier: shuffle conditions around map

%% setup and load data
v = 14;
rfMappingNewMode = 1;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('RF Mapping Analysis - LFP\n');
fprintf('RF Mapping Mode: %d\n', rfMappingNewMode);
fprintf('Session index: %d\n', sessionInd);
fprintf('LFP Channels to Load: %d\n', channelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('RF Mapping info file name: %s\n', rfMappingNewInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

%% input check
assert(numel(channelsToLoad) == 32);

%% load recording information
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, channelsToLoad, 'RFM_NEW', 'LFP_RFM', 0, 1, 1, 0, rfMappingNewInfoFileName, rfMappingNewMode);
sessionName = R.sessionName;
areaName = R.areaName;

plotFileNamePrefix = sprintf('%s-ind%d-%s-ch%d-ch%d-%s', sessionName, sessionInd, areaName, channelsToLoad([1 end]), blockName);
    
%%
trialStructs = loadjson(sprintf('%s/%s', R.rfmResultsRootDir, R.rfmResultsFileNames{1})); % for now, just the first file

clear flashParams;
flashCount = 0;
for i = 1:numel(trialStructs)
    for j = 2:numel(trialStructs{i}.events)
        flashCount = flashCount + 1;
        if flashCount == 1
            flashParams = trialStructs{i}.events{j}.flashParams;
        else
            flashParams(flashCount) = trialStructs{i}.events{j}.flashParams;
        end
    end
end

%         stimId: 34
%        distToFix: 180
%       polarAngle: 0.4490
%         diameter: 59
%     gratingAngle: 45
%         contrast: 0.5000
%        driftRate: 0

% flashStats(i,1:2) = distance(i) * [cos(theta(i)) sin(theta(i))];
% flashStats(i,3) = gratingAngles(gratingAngleBinaryCode + 1);
% flashStats(i,4:6) = [distBinaryCode + 1 polarAngleBinaryCode gratingAngleBinaryCode + 1];

%%
simultSignalTol = 0.001; % tolerance level for simultaneous signals --
% if two events are within this number of  seconds from each other, they 
% represent the same signal. 
flashDuration = 0.1;
flashDurationTol = 0.025;

allEventTimes = D.events;
flashEventTimes = D.events{6};
flashOnsets = nan(size(flashEventTimes));
stimIDs = nan(size(flashEventTimes));
allEventTimes{15} = -1;
allEventTimes{16} = -1; % TODO grating angle doesn't change. there are no EVT15, EVT16
for i = 1:numel(flashEventTimes)
    isFoundEvt = 0;
    foundEvts = cell(16, 1);
    for j = 9:16
        foundEvts{j} = find(abs(flashEventTimes(i) - allEventTimes{j}) < simultSignalTol);
        % there should be at most 1 event9, event10, etc within the
        % simultaneous signal tolerance of each flashEventTime
        assert(numel(foundEvts{j}) <= 1)
        if ~isempty(foundEvts{j})
            isFoundEvt = 1;
        end
    end
    if ~isFoundEvt
        fprintf('huh... no START event 9-16 (flash info) corresponding to event 6 (flash onset) for flash event #%d\n', i);
        continue;
    end

    allFoundEvtTimes = [allEventTimes{9}(foundEvts{9}); 
            allEventTimes{10}(foundEvts{10}); 
            allEventTimes{11}(foundEvts{11}); 
            allEventTimes{12}(foundEvts{12}); 
            allEventTimes{13}(foundEvts{13});
            allEventTimes{14}(foundEvts{14}); 
            allEventTimes{15}(foundEvts{15}); 
            allEventTimes{16}(foundEvts{16})];
    % TODO: do pairwise comparisons instead of consecutive comparisons
    % but this should be good enough
    assert(all(diff(allFoundEvtTimes) < simultSignalTol));
    flashOnset = allFoundEvtTimes(1);
    
    % for confirmation:
    % there should be a second event (in events 9-16) marking the flash
    % offset, indicating the trial was not aborted mid-flash
    isFoundEndEvt = 0;
    foundEndEvts = cell(16, 1);
    for j = 9:16
        foundEndEvts{j} = find(abs(allEventTimes{j} - flashOnset - flashDuration) < flashDurationTol);
        % there should be at most 1 event9, event10, etc within the
        % simultaneous signal tolerance of each flashEventTime
        if numel(foundEndEvts{j}) > 1
            warning('Found %d offset events for event #%d after flash onset #%d at time %f\n', ...
                    numel(foundEndEvts{j}), j, i, flashEventTimes(i));
            foundEndEvts{j} = foundEndEvts{j}(1); % take the earliest if there are multiple
        end
        if ~isempty(foundEndEvts{j})
            isFoundEndEvt = 1;
        end
    end
    if ~isFoundEndEvt
        fprintf('huh... no END event 9-16 (flash info) corresponding to event 6 (flash onset) for flash event #%d\n', i);
        continue;
    end
    allFoundEndEvtTimes = [allEventTimes{9}(foundEndEvts{9}); 
            allEventTimes{10}(foundEndEvts{10}); 
            allEventTimes{11}(foundEndEvts{11}); 
            allEventTimes{12}(foundEndEvts{12}); 
            allEventTimes{13}(foundEndEvts{13});
            allEventTimes{14}(foundEndEvts{14}); 
            allEventTimes{15}(foundEndEvts{15}); 
            allEventTimes{16}(foundEndEvts{16})];
    % TODO: do pairwise comparisons instead of consecutive comparisons
    % but this should be good enough
    assert(all(diff(allFoundEndEvtTimes) < simultSignalTol));
    flashOffset = allFoundEndEvtTimes(1);
    
    flashOnsets(i) = flashOnset;
    stimIDs(i) = bi2de(cellfun(@isempty, foundEvts(9:16))' == 0);
end

flashesToSkip = isnan(stimIDs);
fprintf('Skipping %d flashes due to missing events.\n', sum(flashesToSkip));

%% quick hack
fpStimIDs = cell2mat({flashParams.stimId});
fpCount = numel(fpStimIDs);

nGoodStimIDs = numel(stimIDs(~flashesToSkip));
if nGoodStimIDs < fpCount
    warning('Extra entry/entries in flashParams (JSON) not found in EVT06');
    stimIDsPlus = nan(fpCount,1);
    stimIDsPlus(1:nGoodStimIDs) = stimIDs(~flashesToSkip);
    [(1:fpCount)' stimIDsPlus fpStimIDs'] % print
    
    % remove extra entries in flashParams
    for i = numel(stimIDs):-1:1
        if stimIDs(i) == fpStimIDs(fpCount)
            fpCount = fpCount - 1;
            continue;
        end
        if stimIDs(i) == fpStimIDs(fpCount - 1)
            flashParams(fpCount) = [];
            fpCount = fpCount - 2;
        end
    end
    fpStimIDs = cell2mat({flashParams.stimId});
    fpCount = numel(fpStimIDs);

elseif nGoodStimIDs > fpCount
    fpStimIDsPlus = nan(nGoodStimIDs, 1);
    fpStimIDsPlus(1:fpCount) = fpStimIDs;
    [(1:nGoodStimIDs)' stimIDs(~flashesToSkip) fpStimIDsPlus] % print
end

fprintf('Found %d flash event times and %d flashParam stimIDs.\n', ...
        numel(flashEventTimes) - sum(flashesToSkip), fpCount);
assert((numel(flashEventTimes) - sum(flashesToSkip)) == fpCount);

flashOnsets(flashesToSkip) = [];
stimIDs(flashesToSkip) = [];

%% preprocess LFPs
Fs = D.lfpFs;
nChannels = D.nLfpCh;

% plotNumSpikesByChannel(channelsToLoad, D.allMUAStructs, processedDataDir, blockName, v);
D.adjLfpsClean = interpolateLfpOverSpikeTimes(D.adjLfps, channelsToLoad, Fs, D.allMUAStructs);

hiCutoffFreq = 10;
[channelDataCARNorm,~,~,isNoisyChannel] = preprocessLfps(...
        D.adjLfpsClean, Fs, D.lfpNames, processedDataDir, plotFileNamePrefix, hiCutoffFreq, 1, v);
D.adjLfpsClean = [];

% channelDataBIPNorm = channelDataCARNorm(2:end,:) - channelDataCARNorm(1:end-1,:);
% note channelDataCARNorm, channelDataNorm, and channelDataBIPNorm are HUGE

outlierCheckWindowOffset = [-0.25 0.3];
[flashEventsClean,isEventOutlier] = detectOutlierLfpEvents(channelDataCARNorm, Fs, D.lfpNames, flashOnsets, ...
        outlierCheckWindowOffset, processedDataDir, plotFileNamePrefix, 1, v);
nFlashes = numel(flashEventsClean);

flashParamsClean = flashParams(~isEventOutlier);
distsToFix = cell2mat({flashParamsClean.distToFix});
polarAngles = cell2mat({flashParamsClean.polarAngle});
diameters = cell2mat({flashParamsClean.diameter});
assert(numel(flashEventsClean) == numel(flashParamsClean));

distsToFixUnique = sort(unique(distsToFix));
polarAnglesUnique = sort(unique(polarAngles));

%% process RF for each channel
periEventWindowOffset = [-0.25 0.25];
analysisWindowOffset = [0.025 0.2];
baselineWindowOffset = [-0.175 0];

% set time vector and logical vector for baseline and analysis windows
t = periEventWindowOffset(1):1/Fs:periEventWindowOffset(2)-1/Fs;
analysisWindowLogical = getTimeLogicalWithTolerance(t, analysisWindowOffset);
baselineWindowLogical = getTimeLogicalWithTolerance(t, baselineWindowOffset);

nTime = numel(t);
nEcc = numel(distsToFixUnique);
nPolarAngle = numel(polarAnglesUnique);

for j = 1:nChannels

channelName = D.lfpNames{j};
fprintf('Processing %s (%d/%d = %d%%)... \n', channelName, j, ...
        nChannels, round(j/nChannels*100));

% initialize matrices
averageFlashResponse = zeros(nEcc, nPolarAngle, nTime);
averageFlashResponseNorm = zeros(nEcc, nPolarAngle, nTime);
meanMeanBaselineResponse = nan(nEcc, nPolarAngle);
sdMeanBaselineResponse = nan(nEcc, nPolarAngle);
flashResponseMetric = zeros(nEcc, nPolarAngle);
baselineResponseMetric = zeros(nEcc, nPolarAngle);
trialCounts = zeros(nEcc, nPolarAngle);

for k = 1:nEcc
    for l = 1:nPolarAngle
        % pull all flash events with this eccentricity and polar 
        % angle, regardless of orientation
        theseFlashOnsets = flashEventsClean(distsToFix == distsToFixUnique(k) & polarAngles == polarAnglesUnique(l));
        trialCounts(k,l) = numel(theseFlashOnsets);
        if isempty(theseFlashOnsets)
            continue;
        end

        rawSignals = nan(numel(theseFlashOnsets), nTime);
        % convert flash time to lfp variable index
        % slightly more accurate than createdatamatc b/c of rounding after offset
        startIndices = round((theseFlashOnsets + periEventWindowOffset(1)) * Fs); 
        endIndices = startIndices + round(diff(periEventWindowOffset) * Fs) - 1;

        for n = 1:numel(theseFlashOnsets)
            rawSignals(n,:) = channelDataCARNorm(j,startIndices(n):endIndices(n));
        end
        averageFlashResponse(k,l,:) = mean(rawSignals, 1); % mean over flashes
        meanMeanBaselineResponse(k,l) = mean(averageFlashResponse(k,l,baselineWindowLogical));
        sdMeanBaselineResponse(k,l) = std(averageFlashResponse(k,l,baselineWindowLogical));
%         fprintf('Location (%d,%0.2f): Bsaeline Mean: %0.2f, Baseline SD: %0.2f\n', distsToFixUnique(k), polarAnglesUnique(l), ...
%                 meanMeanBaselineResponse(k,l), sdMeanBaselineResponse(k,l));
        
        % per-condition z-score normalization of LFP
        averageFlashResponseNorm(k,l,:) = (averageFlashResponse(k,l,:) - meanMeanBaselineResponse(k,l)) / sdMeanBaselineResponse(k,l);
        
        % RESPONSE METRIC: max - min in window
        % alts: mean(abs(x)), power(x) at some freq
        flashResponseMetric(k,l) = max(averageFlashResponseNorm(k,l,analysisWindowLogical)) - ...
                min(averageFlashResponseNorm(k,l,analysisWindowLogical));
        baselineResponseMetric(k,l) = max(averageFlashResponseNorm(k,l,baselineWindowLogical)) - ...
                min(averageFlashResponseNorm(k,l,baselineWindowLogical));
    end
end

%% extract all raw signals regardless of flash parameters
rawSignalsAll = nan(numel(flashEventsClean), nTime);
% convert flash time to lfp variable index
% slightly more accurate than createdatamatc b/c of rounding after offset
startIndices = round((flashEventsClean + periEventWindowOffset(1)) * Fs); 
endIndices = startIndices + round(diff(periEventWindowOffset) * Fs) - 1;

for n = 1:numel(flashEventsClean)
    rawSignalsAll(n,:) = channelDataCARNorm(j,startIndices(n):endIndices(n));
end

%%
% figure TODOs
% add separate colorbars for big and for mini, or use same color axis?
% merge x axis for raw data
% label groups of raw data
% y tick off
f = figure_tr_inch(13, 7.5); clf
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig,'Title'), 'Visible', 'on')

modTitle = sprintf('%s (%d flashes)', channelName, nFlashes);
titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(modTitle, 'FontSize', 14, titleParams{:});

%% location params
exampleLfpsByTimeW = 0.1;
meanLfpsGroupedW = 0.11;
heatW = 0.6;
exampleLfpsByTimeH = 0.5;
meanLfpsGroupedH = 0.85;
heatH = 0.85;

exampleLfpsByTimeLeft1 = 0.05;
meanLfpsGroupedLeft1 = exampleLfpsByTimeLeft1 + exampleLfpsByTimeW + 0.03;
meanLfpsGroupedLeft2 = meanLfpsGroupedLeft1 + meanLfpsGroupedW + 0.02;
heatLeft = meanLfpsGroupedLeft2 + meanLfpsGroupedW;
btm = 0.07;
exampleLfpsByTimeBtm = btm + 0.25;

%% plot example raw data around flash onset for N flashes
axes('Position', [exampleLfpsByTimeLeft1 exampleLfpsByTimeBtm exampleLfpsByTimeW exampleLfpsByTimeH]); 
hold on;

nRawToShow = 20;
indRawToShow = round(linspace(1, size(rawSignalsAll, 1), nRawToShow));
rasterY = 0;
yScale = 0.5;
lineParams = {'Color', [0 0 0], 'LineWidth', 1};
for k = 1:numel(indRawToShow)
    rasterY = rasterY + 1;
    plot(t, rawSignalsAll(indRawToShow(k),:)*yScale + rasterY, lineParams{:});
end

yBounds = [0 rasterY+1];
plot([0 0], yBounds, 'k');
plot([analysisWindowOffset(1) analysisWindowOffset(1)], yBounds, '-', 'Color', [0 0.7 0]);
plot([analysisWindowOffset(2) analysisWindowOffset(2)], yBounds, '-', 'Color', [0 0.7 0]);
xlim([-0.1 0.25]);
ylim(yBounds);
set(gca, 'YDir', 'reverse');
xlabel('Time from flash onset (s)'); 
ylabel('Flash number');
title('Example LFPs by time');

%% info
textParams = {'Units', 'normalized', 'FontSize', 8, 'Interpreter', 'none'};
text(axBig, -0.03, 0.11, {''}, ...
        textParams{:});

%% plot mean evoked potentials organized by eccentricity
if nEcc > 1
    axes('Position', [meanLfpsGroupedLeft1 btm meanLfpsGroupedW meanLfpsGroupedH]); 
    hold on;
    channelSep = 1;
    groupSep = (nPolarAngle - 1) / (nEcc - 1) * 2;
    count1 = 1;
    yScale = 0.25;
    cols = lines(nPolarAngle);
    for k = 1:nEcc
        for l = 1:nPolarAngle
            count1 = count1 + channelSep;
            plot(t, squeeze(averageFlashResponseNorm(k,l,:))*yScale + count1, 'Color', cols(l,:));
        end
        if k < nEcc
            count1 = count1 + groupSep;
        end
    end
    count1 = count1 + 2;
    plot([0 0], [0 count1], '-', 'Color', [0.3 0.3 0.3]);
    plot([analysisWindowOffset(1) analysisWindowOffset(1)], [0 count1], '-', 'Color', [0 0.7 0]);
    plot([analysisWindowOffset(2) analysisWindowOffset(2)], [0 count1], '-', 'Color', [0 0.7 0]);
    xlim([-0.1 0.25]);
    ylim([0 count1]);
    set(gca, 'YTickLabel', {});
    xlabel('Time from flash onset (s)');
    title('VEP by eccentricity');
end

%% plot mean evoked potentials organized by polar angle
if nPolarAngle > 1
    axes('Position', [meanLfpsGroupedLeft2 btm meanLfpsGroupedW meanLfpsGroupedH]); 
    hold on;
    channelSep = 1;
    groupSep = 2;
    count2 = 1;
    yScale = 0.25;
    cols = lines(nPolarAngle);
    for l = 1:nPolarAngle
        for k = 1:nEcc
            count2 = count2 + channelSep;
            plot(t, squeeze(averageFlashResponseNorm(k,l,:))*yScale + count2, 'Color', cols(l,:));
        end
        if l < nPolarAngle
            count2 = count2 + groupSep;
        end
    end
    count2 = count2 + 2;
    plot([0 0], [0 count2], '-', 'Color', [0.3 0.3 0.3]);
    plot([analysisWindowOffset(1) analysisWindowOffset(1)], [0 count2], '-', 'Color', [0 0.7 0]);
    plot([analysisWindowOffset(2) analysisWindowOffset(2)], [0 count2], '-', 'Color', [0 0.7 0]);
    xlim([-0.1 0.25]);
    ylim([0 count2]);
    set(gca, 'YTickLabel', {});
    % xlabel('Time from flash onset (s)');
    title('VEP by polar angle');
    if nEcc > 1
        assert(count1 == count2); % i.e. the y axes are the same for the two raw plots
    end
end

%% make heatmap
mapScale = 1/10;
mapDim = round(800 * mapScale) + 1; % x and y: -400 to 400 
rfmapCount = zeros(mapDim, mapDim);
rfmapSmoothAll = zeros(mapDim, mapDim);
mapXOffset = (mapDim + 1)/2;
mapYOffset = (mapDim + 1)/2;
numPixelsPerDegree = 28;

% mean response across conditions with per condition baseline correction
meanBaselineResponses = mean(baselineResponseMetric(:)); % mean baseline response across conditions
sdBaselineResponses = std(baselineResponseMetric(:)); % sd baseline response across conditions

% baseline correct and divide by overall SD of baseline responses across
% conditions in order to account for channel-wise variability and allow for
% comparisons across channels or sessions
flashResponseMetricNorm = (flashResponseMetric - baselineResponseMetric) / sdBaselineResponses;
meanAllFlashResponsesNorm = mean(flashResponseMetricNorm(:));
sdAllFlashResponsesNorm = std(flashResponseMetricNorm(:));
maxAllFlashResponsesNorm = max(flashResponseMetricNorm(:));

threshold = 0;
thresholdName = 'Response == Baseline';

gaussianFilterExtent = mapDim;
gaussianFiltersByDistToFix = cell(nEcc, 1);
for k = 1:nEcc
    % maintain separate maps for each distsToFix because
    % each distToFix has a unique diameter
    matchingDiameters = diameters(distsToFix == distsToFixUnique(k));
    assert(all(matchingDiameters == matchingDiameters(1)))
    diskDiameter = matchingDiameters(1) * mapScale;
    
%     uniformFilter = fspecial('disk', diskDiameter/2);
%     uniformFilter = uniformFilter/max(uniformFilter(:));
    
    gaussianFilterSigma = diskDiameter / 2;
    gaussianFilter = fspecial('gaussian', [gaussianFilterExtent gaussianFilterExtent], ...
            gaussianFilterSigma);
    gaussianFiltersByDistToFix{k} = gaussianFilter / max(gaussianFilter(:)); % normalize so max is 1
    
     % better way: gaussian filter after uniform filter
    
%     filterExtent = diskDiameter + floor(68 * mapScale); 
%     noisyEdgeUniformFilter = zeros(filterExtent, filterExtent);
%     startUniformFilterInd = round((filterExtent-size(uniformFilter, 1))/2) + 1;
%     uniformFilterInds = startUniformFilterInd:startUniformFilterInd+size(uniformFilter, 1)-1;
%     noisyEdgeUniformFilter(uniformFilterInds,uniformFilterInds) = uniformFilter;
    
%     gaussianFilter = fspecial('gaussian', [34 34]*2, 34);
%     noisyEdgeUniformFilter2 = conv2(, gaussianFilter, 'same');
%     noisyEdgeUniformFilter2(noisyEdgeUniformFilter == 1) = 0;
%     noisyEdgeUniformFilter2 = noisyEdgeUniformFilter2 / max(noisyEdgeUniformFilter2(:));
%     noisyEdgeUniformFilter3 = noisyEdgeUniformFilter + noisyEdgeUniformFilter2;
%     noisyEdgeUniformFilter3(noisyEdgeUniformFilter3 > 1) = 1;
%     
%     figure; imagesc(noisyEdgeUniformFilter3); colorbar;
%     figure; plot(noisyEdgeUniformFilter3(floor(size(noisyEdgeUniformFilter3, 1)/2)+1,:))
end

for k = 1:nEcc
    for l = 1:nPolarAngle
        if trialCounts(k,l) > 0
            flashX = distsToFixUnique(k) * cos(polarAnglesUnique(l));
            flashY = distsToFixUnique(k) * sin(polarAnglesUnique(l));
            mapX = round(flashX * mapScale) + mapXOffset;
            mapY = round(flashY * mapScale) + mapYOffset;
            rfmap = zeros(mapDim, mapDim);
            rfmap(mapX,mapY) = flashResponseMetricNorm(k,l); % can be negative
            rfmapCount(mapX,mapY) = trialCounts(k,l);
            rfmapSmooth = imfilter(rfmap, gaussianFiltersByDistToFix{k}, 'replicate');
            % add to existing smooth map
            rfmapSmoothAll = rfmapSmoothAll + rfmapSmooth;
            clear rfmap rfmapSmooth;
        end
    end
end
peakSmoothMapValue = max(rfmapSmoothAll(:));

%% permutation test, shuffling condition locations in map
nRandomizations = 2;
nullDistPeakSmoothMapValues = nan(nRandomizations, 1);
fprintf('Generating null distribution...\n');
for m = 1:nRandomizations
    randpermk = randperm(nEcc);
    randperml = randperm(nPolarAngle);
    rfmapSmoothAllRand = zeros(mapDim, mapDim);
    for k = 1:nEcc
        for l = 1:nPolarAngle
            if trialCounts(k,l) > 0
                flashX = distsToFixUnique(randpermk(k)) * cos(polarAnglesUnique(randperml(l)));
                flashY = distsToFixUnique(randpermk(k)) * sin(polarAnglesUnique(randperml(l)));
                mapX = round(flashX * mapScale) + mapXOffset;
                mapY = round(flashY * mapScale) + mapYOffset;
                rfmap = zeros(mapDim, mapDim);
                rfmap(mapX,mapY) = flashResponseMetricNorm(k,l); % can be negative
%                 rfmapCount(mapX,mapY) = trialCounts(k,l);
                rfmapSmooth = imfilter(rfmap, gaussianFiltersByDistToFix{randpermk(k)}, 'replicate');
                % add to existing smooth map
                rfmapSmoothAllRand = rfmapSmoothAllRand + rfmapSmooth;
                clear rfmap rfmapSmooth;
            end
        end
    end
    nullDistPeakSmoothMapValues(m) = max(rfmapSmoothAllRand(:));
end
fprintf('\tP-value of seeing peak = %0.2f given null distribution (N=%d): %0.3f\n', ...
        peakSmoothMapValue, nRandomizations, sum(nullDistPeakSmoothMapValues > peakSmoothMapValue) / nRandomizations);

%% plot smoothed heatmap for all grating orientations combined
heatAx = axes('Position', [heatLeft btm heatW heatH]);
hold on;
plotRfMapSmooth(rfmapSmoothAll, rfmapCount, numPixelsPerDegree, mapScale);

title(sprintf('VEP Heatmap (%d-%d ms after flash onset; all orientations)', round(analysisWindowOffset * 1000)), 'Interpreter', 'none');
textParams = {'Units', 'normalized', 'FontSize', 8, 'Color', [1 1 1]};
text(0.02, 0.1, sprintf('Threshold: %s', thresholdName), textParams{:});
text(0.02, 0.08, sprintf('Baseline diff, SD across cond: %0.2f', sdBaselineResponses), textParams{:});
text(0.02, 0.06, sprintf('Response diff norm, mean across cond: %0.2f', meanAllFlashResponsesNorm), textParams{:});
text(0.02, 0.04, sprintf('Response diff norm, SD across cond: %0.2f', sdAllFlashResponsesNorm), textParams{:});
text(0.02, 0.02, sprintf('Max response diff norm: %0.2f', maxAllFlashResponsesNorm), textParams{:});
axis(heatAx, 'square'); % shouldn't do much if dims were set properly
colorbar;

% set heatmap scale
caxis([0 35]);%6*sdAllFlashResponses]); 

%% if channel has outlier SD, make background red
if isNoisyChannel(j)
    set(gcf, 'Color', [1 0.5 0.5]);
end

%% save
plotFileName = sprintf('%s/%s-%s-rfm_mode%d_v%d.png', processedDataDir, channelName, blockName, rfMappingNewMode, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
close;

end