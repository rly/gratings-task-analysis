% function lfpRFMapping(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, channelsToLoad, isEighthsOnly)
% LFP RF Mapping, all channels on a probe
% can't really do one channel at a time because of Common Average
% Referencing
% for a couple of sessions and blocks, I used big stimuli in eighths of 
% the screen to do RF mapping 

%% setup and load data
clear;
v = 12;
sessionName = 'M20170201';
sessionInd = 3;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\';
muaDataDirRoot = ['C:\Users\Ryan\Documents\MATLAB\gratings-task-data\' sessionName];
suaMuaDataDirRoot = muaDataDirRoot;
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
channelsToLoad = 1:32;
muaChannelsToLoad = 1:32;
lfpChannelsToLoad = 1:32;
lfpChannels = lfpChannelsToLoad;
isZeroDistractors = 0;
numRandomizations = 2;
isLoadSortedSua = 0;
isLoadMua = 1;
isEighthsOnly = 0;

channelsToProcess = 26:32;
outputDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\GRC';

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
        sessionInd, channelsToLoad, taskName, 'LFP_RFM', isLoadSortedSua, isLoadMua, 1);
sessionName = R.sessionName;
areaName = R.areaName;

plotFileNamePrefix = sprintf('%s-ind%d-%s-ch%d-ch%d-%s', sessionName, sessionInd, areaName, channelsToLoad([1 end]), blockName);

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

%% preprocess LFPs
Fs = D.lfpFs;
nChannels = D.nLfpCh;

D.adjLfpsClean = interpolateLfpOverSpikeTimes(D.adjLfps, channelsToLoad, Fs, D.allMUAStructs);

hiCutoffFreq = 10;
[channelDataCARNorm,channelDataNorm,commonAverageNorm,isNoisyChannel] = preprocessLfps(...
        D.adjLfpsClean, Fs, D.lfpNames, processedDataDir, plotFileNamePrefix, hiCutoffFreq, 1, v);
D.adjLfpsClean = [];

channelDataBIPNorm = channelDataCARNorm(2:end,:) - channelDataCARNorm(1:end-1,:);
% note channelDataCARNorm, channelDataNorm, and channelDataBIPNorm are HUGE

outlierCheckWindowOffset = [-0.25 0.3];
[flashEventsClean,isEventOutlier] = detectOutlierLfpEvents(channelDataCARNorm, Fs, D.lfpNames, origFlashEvents, ...
        outlierCheckWindowOffset, processedDataDir, plotFileNamePrefix, 1, v);
nFlashes = numel(flashEventsClean);

flashParamsClean = flashParams(~isEventOutlier,:);
assert(numel(flashEventsClean) == size(flashParamsClean, 1));

%% process RF for each channel
thresh = 8;

f = figure_tr_inch(3.5, 3.5); 
clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

heatAx = axes('Position', [0.08 0.08 0.84 0.84]); 
hold on;

cols = hsv(numel(channelsToProcess));

for i = 1:numel(channelsToProcess)
    j = channelsToProcess(i);

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
        theseFlashOnsets = flashEventsClean(flashParamsClean(:,4) == k & flashParamsClean(:,5) == l);
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
            rawSignals(n,:) = channelDataCARNorm(j,startIndices(n):endIndices(n));
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

rawSignals = nan(numel(flashEventsClean), numel(t));
% convert flash time to lfp variable index
% slightly more accurate than createdatamatc b/c of rounding after offset
startIndices = round((flashEventsClean + periEventWindowOffset(1)) * Fs); 
endIndices = startIndices + round(diff(periEventWindowOffset) * Fs) - 1;

for n = 1:numel(flashEventsClean)
    rawSignals(n,:) = channelDataCARNorm(j,startIndices(n):endIndices(n));
end

%% make heatmap
mapScale = 1/2;
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
%     for m = 1:numel(gratingAnglesUnique) + 1
        rfmapSmoothByOri(k,m2,:,:) = imfilter(squeeze(rfmapByOri(k,m2,:,:)), gaussianFilter, 'replicate');
%     end
end

% rfmapAllOri = squeeze(sum(rfmapByOri(:,m2,:,:), 1));
rfmapAllOriCount = squeeze(sum(rfmapByOriCount(:,m2,:,:), 1));

% rfmapSmoothByOriAll = squeeze(sum(rfmapSmoothByOri, 1)); % sum over the different smooth maps for each distsToFix
% rfmapByOriCountAll = squeeze(sum(rfmapByOriCount, 1)); % sum over the different smooth maps for each distsToFix
rfmapAllSmoothByOriAll = squeeze(sum(rfmapSmoothByOri(:,m2,:,:), 1));

%% plot smoothed heatmap for all and individual orientations

% start with all orientations together
hold on;
% plotRfMapSmooth(rfmapAllSmoothByOriAll, rfmapAllOriCount, numPixelsPerDegree, mapScale, mapDim, mapXOffset, mapYOffset);
threshGroups = bwlabel(rfmapAllSmoothByOriAll > thresh, 4);
threshBoundaries = bwboundaries(threshGroups);
x = ((1:mapDim) - mapXOffset) / numPixelsPerDegree / mapScale;
y = ((1:mapDim) - mapYOffset) / numPixelsPerDegree / mapScale;
for k = 1:length(threshBoundaries)
   boundary = threshBoundaries{k};
   plot(x(boundary(:,1)), y(boundary(:,2)), 'Color', cols(i,:), 'LineWidth', 2)
end

xlim(([1 mapDim] - mapXOffset) / numPixelsPerDegree / mapScale);
ylim(([1 mapDim] - mapYOffset) / numPixelsPerDegree / mapScale);

set(gca, 'XTickLabel', []);
set(gca, 'yTickLabel', []);
set(gca, 'box', 'on');
set(gca, 'LineWidth', 1);
set(gca, 'Color', 0.5*ones(3, 1));
axis(gca, 'square')
drawnow;

end

% draw sampled points
[rfmapX, rfmapY] = find(rfmapAllOriCount > 0);
h = plot((rfmapX - mapXOffset) / numPixelsPerDegree / mapScale, ...
        (rfmapY - mapYOffset) / numPixelsPerDegree / mapScale, ...
        '.', 'Color', [0 0 0], 'MarkerSize', 10);
uistack(h, 'bottom');

% draw fixation point
fixW = 1.8;
fixH = fixW;
rectangle('Position', [-fixW/2 -fixH/2 fixW fixH], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);


%% save
plotFileName = sprintf('%s/%s-ch%d-ch%d-%s-v%d.png', outputDir, sessionName, channelsToProcess([1 end]), blockName, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');