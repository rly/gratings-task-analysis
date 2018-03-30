% function lfpRFMapping(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, channelsToLoad, isEighthsOnly)
% LFP RF Mapping, all channels on a probe
% can't really do one channel at a time because of Common Average
% Referencing
% for a couple of sessions and blocks, I used big stimuli in eighths of 
% the screen to do RF mapping 

%% setup and load data
v = 11;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('RF Mapping Analysis - LFP\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('LFP Channel to Load: %d\n', channelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

%% input check
% assert(numel(channelsToLoad) == 1);
assert(numel(channelsToLoad) > 1);

%% load recording information
if ~isEighthsOnly
    taskName = 'RFM_OLD';
else
    taskName = 'RFM_EIGHTHS';
end
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, channelsToLoad, taskName, 'LFP_RFM', 1, 1);

%%
% task as of 1/10/17 (or earlier)
% 325ms fixation
% 25ms of pre-flashes period (continued fixation)
% up to N=4 flashes, each 100ms duration, separated by 200ms ISI
% 250-500ms after 5th flash until FP dims
% 260-775ms after FP dims is correct response window
% 4 distances to fixation 150:60:330
% 11 polar angles (-5:5)*pi/8
% 4 grating orientations (0:3)*pi/4
% 3 reps each

% EVT05 - begin pre-cue period
% EVT06 - flash onset
% EVT07 - FP dim
% EVT08 - juice onset

origFlashEvents = D.events{6};

% update based on the PCL code
if ~isEighthsOnly
    % 4 possible, coded in events 9-10: 01, 01, 10, 11, in the order below
    distsToFixUnique = 150:60:330;
    % 7 possible, coded in events 11-14: 0001, 0010, ..., 1011, in the order below
    polarAnglesUnique = (-5:5)*pi/8;
else
    distsToFixUnique = 240;
    polarAnglesUnique = (-6:2:8)*pi/8;
end

% 4 possible, coded in events 15-16: 00, 01, 10, 11, in the order below
gratingAnglesUnique = (0:3)*pi/4;

[flashOnsets,flashParams] = decodeFlashParams(origFlashEvents, ...
        D.events, distsToFixUnique, polarAnglesUnique, gratingAnglesUnique);
nFlashes = numel(flashOnsets);

%% preprocess LFPs
Fs = D.lfpFs;
nChannels = D.nLfpCh;

D.adjLfpsClean = interpolateLfpOverSpikeTimes(D.adjLfps, channelsToLoad, Fs, D.allMUAStructs);

hiCutoffFreq = 50;
outlierCheckWindowOffset = [-0.25 0.3];
[channelDataNorm,flashOnsetsClean,isEventOutlier,isNoisyChannel] = preprocessLfps(D.adjLfpsClean, ...
        Fs, D.lfpNames, flashOnsets, processedDataDir, blockName, hiCutoffFreq, 1, v, outlierCheckWindowOffset);
D.adjLfps = [];
D.adjLfpsClean = [];
for j = 1:nChannels
    fprintf('Channel %d: Mean: %0.2f, SD: %0.2f\n', j, nanmean(channelDataNorm(j,:)), nanstd(channelDataNorm(j,:)));
end

flashParamsClean = flashParams(~isEventOutlier,:);
assert(numel(flashOnsetsClean) == size(flashParamsClean, 1));

%% process RF for each channel
for j = 1:nChannels

channelName = D.lfpNames{j};
fprintf('Processing %s (%d/%d = %d%%)... \n', channelName, j, ...
        nChannels, round(j/nChannels*100));

periEventWindowOffset = [-0.25 0.25];
t = periEventWindowOffset(1):1/Fs:periEventWindowOffset(2)-1/Fs;

analysisWindowOffset = [0.025 0.2];
analysisWindowLogical = t >= analysisWindowOffset(1) & t < analysisWindowOffset(2);

baselineWindowOffset = [-0.175 0];
baselineWindowLogical = t >= baselineWindowOffset(1) & t < baselineWindowOffset(2);

averageFlashResponse = zeros(numel(distsToFixUnique), numel(polarAnglesUnique), numel(gratingAnglesUnique)+1, numel(t));
averageFlashResponseNorm = zeros(numel(distsToFixUnique), numel(polarAnglesUnique), numel(gratingAnglesUnique)+1, numel(t));
flashResponseMetric = zeros(numel(distsToFixUnique), numel(polarAnglesUnique), numel(gratingAnglesUnique)+1);
baselineResponseMetric = zeros(numel(distsToFixUnique), numel(polarAnglesUnique), numel(gratingAnglesUnique)+1);
trialCounts = zeros(numel(distsToFixUnique), numel(polarAnglesUnique), numel(gratingAnglesUnique)+1);
meanMeanBaselineResponse = zeros(numel(distsToFixUnique), numel(polarAnglesUnique), numel(gratingAnglesUnique)+1);
sdMeanBaselineResponse = zeros(numel(distsToFixUnique), numel(polarAnglesUnique), numel(gratingAnglesUnique)+1);
m2 = numel(gratingAnglesUnique) + 1; % for indexing the no orientation case

% TEMP: use only index m2 into the gratingAngles dimension (3rd dimension
% typically). currently ignoring grating orientation

for k = 1:numel(distsToFixUnique)
    for l = 1:numel(polarAnglesUnique)
%         for m = 1:numel(gratingAnglesUnique)
%             theseFlashOnsets = flashOnsetsClean(flashStatsClean(:,4) == k & flashStatsClean(:,5) == l & flashStatsClean(:,6) == m);
%             trialCounts(k,l,m) = numel(theseFlashOnsets);
%             if isempty(theseFlashOnsets)
%                 continue;
%             end
%                         
%             rawSignals = nan(numel(theseFlashOnsets), numel(t));
%             % convert flash time to lfp variable index
%             % slightly more accurate than createdatamatc b/c of rounding after offset
%             startIndices = round((theseFlashOnsets + periEventWindowOffset(1)) * Fs); 
%             endIndices = startIndices + round(diff(periEventWindowOffset) * Fs) - 1;
%             
%             for n = 1:numel(theseFlashOnsets)
%                 rawSignals(n,:) = channelDataNorm(j,startIndices(n):endIndices(n));
%             end
%             averageFlashResponse(k,l,m,:) = mean(rawSignals, 1); % mean over flashes
% %             averageFlashResponseMeanOverTime(k,l,m) = mean(abs(averageFlashResponse(k,l,m,:)));
%             responseMetric(k,l,m) = max(averageFlashResponse(k,l,m,analysisWindowLogical)) - ...
%                     min(averageFlashResponse(k,l,m,analysisWindowLogical));
%             baselineMetric(k,l,m) = max(averageFlashResponse(k,l,m,baselineWindowLogical)) - ...
%                     min(averageFlashResponse(k,l,m,baselineWindowLogical));
%         end
        
        % case: regardless of orientation 
        theseFlashOnsets = flashOnsetsClean(flashParamsClean(:,4) == k & flashParamsClean(:,5) == l);
        trialCounts(k,l,m2) = numel(theseFlashOnsets);
        if isempty(theseFlashOnsets)
            continue;
        end
        
        rawSignals = nan(numel(theseFlashOnsets), numel(t));
        % convert flash time to lfp variable index
        % slightly more accurate than createdatamatc b/c of rounding after offset
        startIndices = round((theseFlashOnsets + periEventWindowOffset(1)) * Fs); 
        endIndices = startIndices + round(diff(periEventWindowOffset) * Fs) - 1;
        
        for n = 1:numel(theseFlashOnsets)
            rawSignals(n,:) = channelDataNorm(j,startIndices(n):endIndices(n));
        end
        
        averageFlashResponse(k,l,m2,:) = mean(rawSignals, 1); % mean over flashes
        meanMeanBaselineResponse(k,l,m2) = mean(averageFlashResponse(k,l,m2,baselineWindowLogical));
        sdMeanBaselineResponse(k,l,m2) = std(averageFlashResponse(k,l,m2,baselineWindowLogical));
%         fprintf('Location (%d,%0.2f): Bsaeline Mean: %0.2f, Baseline SD: %0.2f\n', distsToFixUnique(k), polarAnglesUnique(l), ...
%                 meanMeanBaselineResponse(k,l), sdMeanBaselineResponse(k,l));
        
        % per-condition z-score normalization of LFP
        averageFlashResponseNorm(k,l,m2,:) = (averageFlashResponse(k,l,m2,:) - meanMeanBaselineResponse(k,l,m2)) / sdMeanBaselineResponse(k,l,m2);
        
%         averageFlashResponseMeanOverTime(k,l,m2) = mean(abs(averageFlashResponse(k,l,m2,:)));
        flashResponseMetric(k,l,m2) = max(averageFlashResponseNorm(k,l,m2,analysisWindowLogical)) - ...
                min(averageFlashResponseNorm(k,l,m2,analysisWindowLogical));
        baselineResponseMetric(k,l,m2) = max(averageFlashResponseNorm(k,l,m2,baselineWindowLogical)) - ...
                    min(averageFlashResponseNorm(k,l,m2,baselineWindowLogical));
    end
end

rawSignals = nan(numel(flashOnsetsClean), numel(t));
% convert flash time to lfp variable index
% slightly more accurate than createdatamatc b/c of rounding after offset
startIndices = round((flashOnsetsClean + periEventWindowOffset(1)) * Fs); 
endIndices = startIndices + round(diff(periEventWindowOffset) * Fs) - 1;

for n = 1:numel(flashOnsetsClean)
    rawSignals(n,:) = channelDataNorm(j,startIndices(n):endIndices(n));
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
waveformW = 0.1;
rasterByTimeW = 0.1;
rawW = 0.11;
heatW = 0.6;
% miniHeatW = 0.18;
waveformH = 0.15;
rasterByTimeH = 0.35;
rawH = 0.85;
heatH = 0.85;
% miniHeatH = 0.385;

waveformLeft1 = 0.05;
rasterByTimeLeft1 = waveformLeft1;
rawLeft1 = rasterByTimeLeft1 + rasterByTimeW + 0.03;
rawLeft2 = rawLeft1 + rawW + 0.02;
heatLeft = rawLeft2 + rawW;
% miniHeatLeft1 = heatLeft + heatW + 0.005;
% miniHeatLeft2 = miniHeatLeft1 + miniHeatW + 0.01;
btm = 0.07;
rasterByTimeBtm = btm + 0.25;
waveformBtm = rasterByTimeBtm + rasterByTimeH + 0.1;
% miniHeatBtm2 = btm + miniHeatH + 0.08;

%%
axes('Position', [rasterByTimeLeft1 rasterByTimeBtm rasterByTimeW rasterByTimeH]); 
hold on;

nRawToShow = 20;
indRawToShow = round(linspace(1, size(rawSignals, 1), nRawToShow));
rasterY = 0;
yScale = 0.5;
lineParams = {'Color', [0 0 0], 'LineWidth', 1};
for k = 1:numel(indRawToShow)
    rasterY = rasterY + 1;
    plot(t, rawSignals(indRawToShow(k),:)*yScale + rasterY, lineParams{:});
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
titleH = title('Example LFPs by time', 'Interpreter', 'None');
set(gca, 'Layer', 'top');


%% info
textParams = {'Units', 'normalized', 'FontSize', 8, 'Interpreter', 'none'};
text(axBig, -0.03, 0.11, {''}, ...
        textParams{:});

%%
if numel(distsToFixUnique) > 1
    axes('Position', [rawLeft1 btm rawW rawH]); 
    hold on;
    channelSep = 1;
    groupSep = (numel(polarAnglesUnique) - 1) / (numel(distsToFixUnique) - 1) * 2;
    count1 = 1;
    yScale = 1;
    cols = lines(numel(polarAnglesUnique));
    for k = 1:numel(distsToFixUnique)
        for l = 1:numel(polarAnglesUnique)
            count1 = count1 + channelSep;
            plot(t, squeeze(averageFlashResponse(k,l,m2,:))*yScale + count1, 'Color', cols(l,:));
        end
        if numel(polarAnglesUnique) > 1 && k < numel(distsToFixUnique)
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

%%
if numel(polarAnglesUnique) > 1
    axes('Position', [rawLeft2 btm rawW rawH]); 
    hold on;
    channelSep = 1;
    groupSep = 2;
    count2 = 1;
    yScale = 1;
    cols = lines(numel(polarAnglesUnique));
    for l = 1:numel(polarAnglesUnique)
        for k = 1:numel(distsToFixUnique)
            count2 = count2 + channelSep;
            plot(t, squeeze(averageFlashResponse(k,l,m2,:))*yScale + count2, 'Color', cols(l,:));
        end
        if numel(distsToFixUnique) > 1 && l < numel(polarAnglesUnique)
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
    if numel(distsToFixUnique) > 1
        assert(count1 == count2); % i.e. the y axes are the same for the two raw plots
    end
end

%% make heatmap
mapScale = 1/10;
mapDim = round(800 * mapScale) + 1; % x and y: -400 to 400 
rfmapByOri = zeros(numel(distsToFixUnique), numel(gratingAnglesUnique)+1, mapDim, mapDim); % extra dim uses all trials
rfmapByOriCount = zeros(numel(distsToFixUnique), numel(gratingAnglesUnique)+1, mapDim, mapDim);
rfmapSmoothByOri = zeros(size(rfmapByOri));
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
for k = 1:numel(distsToFixUnique)
    % slight adjustment for eve movement, which is radius 68 but here let's treat it as diameter 68
    if ~isEighthsOnly
        diskDiameter = floor((distsToFixUnique(k) / 150 * 3 * 28 * 2 - 68) * mapScale); % pixels, floor to nearest even
    else
        diskDiameter = floor((distsToFixUnique(k) * 0.04 * 28 * 2 - 68) * mapScale); % pixels, floor to nearest even
    end
%     uniformFilter = fspecial('disk', diskDiameter/2);
%     uniformFilter = uniformFilter/max(uniformFilter(:));
    
    gaussianFilterSigma = diskDiameter / 3;
    gaussianFilter = fspecial('gaussian', [gaussianFilterExtent gaussianFilterExtent], ...
            gaussianFilterSigma);
    gaussianFilter = gaussianFilter / max(gaussianFilter(:));
    
    for l = 1:numel(polarAnglesUnique)
        flashX = distsToFixUnique(k) * cos(polarAnglesUnique(l));
        flashY = distsToFixUnique(k) * sin(polarAnglesUnique(l));
        mapX = round(flashX * mapScale) + mapXOffset;
        mapY = round(flashY * mapScale) + mapYOffset;
%         for m = 1:numel(gratingAnglesUnique)
%             % threshold at baseline
%             rfmapByOri(k,m,mapX,mapY) = max(0, responseMetric(k,l,m) - threshold);
%             rfmapByOriCount(k,m,mapX,mapY) = trialCounts(k,l,m); 
%         end
%         if sum(trialCounts(k,l,:)) > 0
%             % weighted mean across orientations
%             % threshold at baseline
%             rfmapByOri(k,m2,mapX,mapY) = max(0, sum(trialCounts(k,l,:).*responseMetric(k,l,:)) ...
%                     / sum(trialCounts(k,l,:)) - threshold);
%             rfmapByOriCount(k,m2,mapX,mapY) = sum(trialCounts(k,l,:)); 
%         end
%         rfmapByOri(k,m2,mapX,mapY) = max(0, responseMetric(k,l,m2) - threshold);
        rfmapByOri(k,m2,mapX,mapY) = flashResponseMetricNorm(k,l,m2); % can be negative
        rfmapByOriCount(k,m2,mapX,mapY) = trialCounts(k,l,m2); 
    end
    % maintain separate smooth maps for each distsToFix
    for m = 1:numel(gratingAnglesUnique) + 1
        rfmapSmoothByOri(k,m,:,:) = imfilter(squeeze(rfmapByOri(k,m,:,:)), gaussianFilter, 'replicate');
    end
end

% rfmapAllOri = squeeze(sum(rfmapByOri(:,m2,:,:), 1));
rfmapAllOriCount = squeeze(sum(rfmapByOriCount(:,m2,:,:), 1));

% rfmapSmoothByOriAll = squeeze(sum(rfmapSmoothByOri, 1)); % sum over the different smooth maps for each distsToFix
% rfmapByOriCountAll = squeeze(sum(rfmapByOriCount, 1)); % sum over the different smooth maps for each distsToFix
rfmapAllSmoothByOriAll = squeeze(sum(rfmapSmoothByOri(:,m2,:,:), 1));

%%
% flatFilterSize = 4;
% flatFilter = ones(flatFilterSize, flatFilterSize);
% figure;
% hold on;
% rfmapSmooth = conv2(rfmapAllOri, flatFilter, 'same');
% imagesc(((1:mapDim) - mapXOffset) / numPixelsPerDegree, ...
%         ((1:mapDim) - mapYOffset) / numPixelsPerDegree, rfmapSmooth');
% xlim(([1 mapDim] - mapXOffset) / numPixelsPerDegree);
% ylim(([1 mapDim] - mapYOffset) / numPixelsPerDegree);
% xlabel('Degrees horizontal');
% ylabel('Degrees vertical');
% colorbar;
% 
% title(sprintf('%s - Response Heatmap (20-220ms after flash onset)', spkVarName), 'Interpreter', 'none');

%% plot smoothed heatmap for all and individual orientations

% start with all orientations together
heatAx = axes('Position', [heatLeft btm heatW heatH]); 
hold on;
plotRfMapSmooth(rfmapAllSmoothByOriAll, rfmapAllOriCount, numPixelsPerDegree, mapScale, mapDim, mapXOffset, mapYOffset);
title(sprintf('VEP Heatmap (%d-%d ms after flash onset; all orientations)', round(analysisWindowOffset * 1000)), 'Interpreter', 'none');
textParams = {'Units', 'normalized', 'FontSize', 8, 'Color', [1 1 1]};
text(0.02, 0.1, sprintf('Threshold: %s', thresholdName), textParams{:});
text(0.02, 0.08, sprintf('Baseline diff, SD across cond: %0.2f', sdBaselineResponses), textParams{:});
text(0.02, 0.06, sprintf('Response diff norm, mean across cond: %0.2f', meanAllFlashResponsesNorm), textParams{:});
text(0.02, 0.04, sprintf('Response diff norm, SD across cond: %0.2f', sdAllFlashResponsesNorm), textParams{:});
text(0.02, 0.02, sprintf('Max response diff norm: %0.2f', maxAllFlashResponsesNorm), textParams{:});
axis(heatAx, 'square'); % shouldn't do much if i set the dims properly
colorbar;
caxis([0 15]); 

%%
if isNoisyChannel(j)
    set(gcf, 'Color', [1 0.5 0.5]); % make background red
end

plotFileName = sprintf('%s/%s-%s-v%d.png', processedDataDir, channelName, blockName, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
close;

end