function muaRFMapping(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad)
% MUA RF Mapping, one channel

%% setup and load data
v = 11;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('RF Mapping Analysis - MUA\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('MUA Channel to Load: %d\n', muaChannelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

%% input check
assert(numel(muaChannelsToLoad) == 1);

%% load recording information
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, muaChannelsToLoad, 'RFM_OLD', 'MUA_RFM', 0, 1, 0, 0);

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
nFlashes = numel(origFlashEvents);

% update based on the PCL code
% 4 possible, coded in events 9-10: 01, 01, 10, 11, in the order below
distsToFix = 150:60:330;
% 7 possible, coded in events 11-14: 0001, 0010, ..., 1011, in the order below
polarAngles = (-5:5)*pi/8;
% 4 possible, coded in events 15-16: 00, 01, 10, 11, in the order below
gratingAngles = (0:3)*pi/4;

[flashOnsets,flashStats] = decodeFlashParams(origFlashEvents, ...
        D.events, distsToFix, polarAngles, gratingAngles);

%%
nUnits = numel(D.allMUAStructs);
assert(nUnits == 1);
j = 1;

%%
fprintf('Processing %d MUAs...\n', nUnits);

spikeStruct = D.allMUAStructs{j};
unitName = spikeStruct.name;
spikeTimes = spikeStruct.ts;
fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
        nUnits, round(j/nUnits*100));

%% compute psth for each flash condition
% this does not subtract out per-condition baseline activity
psthWindow = [0.25 0.25]; % secs before, secs after
analysisWindowOffset = [0.02 0.22]; % secs offset from event onset
analysisWindow = analysisWindowOffset + psthWindow(1);
analysisWindowDur = diff(analysisWindowOffset);

baselineWindowOffset = [-0.2 0];
baselineWindow = baselineWindowOffset + psthWindow(1);
baselineWindowDur = diff(baselineWindowOffset);

kernelSigma = 0.01;
nTime = fix(5*sum(psthWindow)/kernelSigma) + 1; % number of time steps. 
psthT = linspace(0, sum(psthWindow), nTime);
t = psthT - psthWindow(1);

allAlignedSpikeTimes = createdatamatpt(spikeTimes, flashOnsets, psthWindow);

% TODO consider removing quiet trials

m2 = numel(gratingAngles) + 1;
allPsths = nan(numel(distsToFix), numel(polarAngles), numel(gratingAngles) + 1, nTime);
trialCounts = zeros(numel(distsToFix), numel(polarAngles), numel(gratingAngles));
responseMetric = zeros(numel(distsToFix), numel(polarAngles), numel(gratingAngles));
baselineMetric = zeros(numel(distsToFix), numel(polarAngles), numel(gratingAngles));
for k = 1:numel(distsToFix)
    for l = 1:numel(polarAngles)
        for m = 1:numel(gratingAngles)
            theseFlashOnsets = flashOnsets(flashStats(:,4) == k & flashStats(:,5) == l & flashStats(:,6) == m);
            trialCounts(k,l,m) = numel(theseFlashOnsets);
            if isempty(theseFlashOnsets)
                continue;
            end

            dataStruct = createdatamatpt(spikeTimes, theseFlashOnsets, psthWindow);

            theseNSpikesPostFlash = zeros(numel(dataStruct), 1);
            theseNSpikesBaseline = zeros(numel(dataStruct), 1);
            % maybe could just collapse the struct and then do the time check
            for n = 1:numel(dataStruct)
                theseNSpikesPostFlash(n) = sum(analysisWindow(1) <= dataStruct(n).times & ...
                    dataStruct(n).times <= analysisWindow(2));
                theseNSpikesBaseline(n) = sum(baselineWindow(1) <= dataStruct(n).times & ...
                    dataStruct(n).times <= baselineWindow(2));
            end
            responseMetric(k,l,m) = sum(theseNSpikesPostFlash) / numel(dataStruct); % average num spikes per trial
            baselineMetric(k,l,m) = sum(theseNSpikesBaseline) / numel(dataStruct); % average num spikes per trial

            R = fixedPsth(dataStruct, kernelSigma, 2, psthT);
            if ~isempty(R)
                allPsths(k,l,m,:) = R;
            end
        end
        
        % case: regardless of orientation 
        theseFlashOnsets = flashOnsets(flashStats(:,4) == k & flashStats(:,5) == l);
        if isempty(theseFlashOnsets)
            continue;
        end

        dataStruct = createdatamatpt(spikeTimes, theseFlashOnsets, psthWindow);
        R = fixedPsth(dataStruct, kernelSigma, 2, psthT);
        if ~isempty(R)
            allPsths(k,l,m2,:) = R;
        end
    end
end

%%
% figure TODOs
% add separate colorbars for big and for mini, or use same color axis?
% merge x axis for raw data
% label groups of raw data
% y tick off
f = figure_tr_inch(13, 7.5); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on')

modTitle = sprintf('Spike RF Mapping: %s (%d flashes)', unitName, nFlashes);
titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(modTitle, 'FontSize', 14, titleParams{:});

%% location params
waveformW = 0.1;
rasterByTimeW = 0.1;
rawW = 0.11;
heatW = 0.6;
waveformH = 0.15;
rasterByTimeH = 0.32;
rawH = 0.84;
heatH = 0.84;

waveformLeft1 = 0.05;
rasterByTimeLeft1 = waveformLeft1;
rawLeft1 = rasterByTimeLeft1 + rasterByTimeW + 0.03;
rawLeft2 = rawLeft1 + rawW + 0.02;
heatLeft = rawLeft2 + rawW + 0.01;
btm = 0.07;
infoTextTop = btm + 0.18;
rasterByTimeBtm = infoTextTop + 0.09;
waveformBtm = 0.76;%rasterByTimeBtm + rasterByTimeH + 0.1;

%% plot spike waveform
axes('Position', [waveformLeft1 waveformBtm waveformW waveformH]); 
plotSpikeWaveform(D, j, 1);

%% plot raster
axes('Position', [rasterByTimeLeft1 rasterByTimeBtm rasterByTimeW rasterByTimeH]); 
hold on;

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
            plot([adjTimes(l) adjTimes(l)], rasterY + [-lineHeight lineHeight]/2, lineParams{:});
        end
    end
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
titleH = title('Raster by time', 'Interpreter', 'None');
set(gca, 'Layer', 'top');

%% info
writeUnitInfo(spikeStruct, axBig, -0.03, infoTextTop);

%% plot VES by eccentricity
axes('Position', [rawLeft1 btm rawW rawH]); 
hold on;
channelSep = 1;
groupSep = (numel(polarAngles) - 1) / (numel(distsToFix) - 1) * 2;
count1 = 1;
yScale = 0.07;
cols = lines(numel(polarAngles));
for k = 1:numel(distsToFix)
    for l = 1:numel(polarAngles)
        count1 = count1 + channelSep;
        plot(t, squeeze(allPsths(k,l,m2,:))*yScale + count1, 'Color', cols(l,:));
    end
    if k < numel(distsToFix)
        count1 = count1 + groupSep;
    end
end
count1 = count1 + 6;
plot([0 0], [0 count1], '-', 'Color', [0.3 0.3 0.3]);
fillH = jbfill(analysisWindowOffset, [-1 -1], [count1 count1], [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
xlim([-0.1 0.25]);
ylim([0 count1]);
set(gca, 'YTickLabel', {});
xlabel('Time from flash onset (s)');
title('VES by eccentricity');

%% plot VES by polar angle
axes('Position', [rawLeft2 btm rawW rawH]); 
hold on;
channelSep = 1;
groupSep = 2;
count2 = 1;
yScale = 0.07;
cols = lines(numel(polarAngles));
for l = 1:numel(polarAngles)
    for k = 1:numel(distsToFix)
        count2 = count2 + channelSep;
        plot(t, squeeze(allPsths(k,l,m2,:))*yScale + count2, 'Color', cols(l,:));
    end
    if l < numel(polarAngles)
        count2 = count2 + groupSep;
    end
end
count2 = count2 + 6;
plot([0 0], [0 count2], '-', 'Color', [0.3 0.3 0.3]);
fillH = jbfill(analysisWindowOffset, [-1 -1], [count2 count2], [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
xlim([-0.1 0.25]);
ylim([0 count2]);
set(gca, 'YTickLabel', {});
% xlabel('Time from flash onset (s)');
title('VES by polar angle');
assert(count1 == count2); % i.e. the y axes are the same for the two raw plots

%% make heatmap of average number of spikes (doesn't use psth)
mapScale = 1/10;
mapDim = round(800 * mapScale) + 1; % x and y: -400 to 400 
rfmapByOri = zeros(numel(distsToFix), numel(gratingAngles)+1, mapDim, mapDim); % extra dim uses all trials
rfmapByOriCount = zeros(numel(distsToFix), numel(gratingAngles)+1, mapDim, mapDim);
rfmapSmoothByOri = zeros(size(rfmapByOri));
mapXOffset = (mapDim + 1)/2;
mapYOffset = (mapDim + 1)/2;
numPixelsPerDegree = 28;

meanAllResponses = mean(responseMetric(:));
baselineResponse = mean(baselineMetric(:));
threshold = baselineResponse;
thresholdName = 'baseline';

gaussianFilterExtent = mapDim;
for k = 1:numel(distsToFix)
    % slight adjustment for eve movement, which is radius 68 but here let's treat it as diameter 68
    diskDiameter = floor((distsToFix(k) / 150 * 3 * 28 * 2 - 68) * mapScale); % pixels, floor to nearest even
    uniformFilter = fspecial('disk', diskDiameter/2);
    uniformFilter = uniformFilter/max(uniformFilter(:));
    
    gaussianFilterSigma = diskDiameter / 3;
    gaussianFilter = fspecial('gaussian', [gaussianFilterExtent gaussianFilterExtent], ...
            gaussianFilterSigma);
    gaussianFilter = gaussianFilter / max(gaussianFilter(:));
    
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
    
    for l = 1:numel(polarAngles)
        flashX = distsToFix(k) * cos(polarAngles(l));
        flashY = distsToFix(k) * sin(polarAngles(l));
        mapX = round(flashX * mapScale) + mapXOffset;
        mapY = round(flashY * mapScale) + mapYOffset;
        for m = 1:numel(gratingAngles)
            % threshold at baseline
            rfmapByOri(k,m,mapX,mapY) = max(0, responseMetric(k,l,m) - threshold);
            rfmapByOriCount(k,m,mapX,mapY) = trialCounts(k,l,m); 
        end
        if sum(trialCounts(k,l,:)) > 0
            % weighted mean across orientations
            % threshold at baseline
            rfmapByOri(k,m2,mapX,mapY) = max(0, sum(trialCounts(k,l,:).*responseMetric(k,l,:)) ...
                    / sum(trialCounts(k,l,:)) - threshold);
            rfmapByOriCount(k,m2,mapX,mapY) = sum(trialCounts(k,l,:)); 
        end
    end
    % maintain separate smooth maps for each distsToFix
    for m = 1:numel(gratingAngles) + 1
        rfmapSmoothByOri(k,m,:,:) = imfilter(squeeze(rfmapByOri(k,m,:,:)), gaussianFilter, 'replicate');
    end
end

rfmapAllOri = squeeze(sum(rfmapByOri(:,numel(gratingAngles)+1,:,:), 1));
rfmapAllOriCount = squeeze(sum(rfmapByOriCount(:,numel(gratingAngles)+1,:,:), 1));

rfmapSmoothByOriAll = squeeze(sum(rfmapSmoothByOri, 1)); % sum over the different smooth maps for each distsToFix
rfmapByOriCountAll = squeeze(sum(rfmapByOriCount, 1)); % sum over the different smooth maps for each distsToFix
rfmapAllSmoothByOriAll = squeeze(sum(rfmapSmoothByOri(:,numel(gratingAngles)+1,:,:), 1));

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

% compute spike rate stats, not using psth
avgSpikeRatePostFlash = responseMetric / analysisWindowDur;
avgSpikeRateBaseline = baselineMetric / baselineWindowDur;
meanSpikeRate = mean(avgSpikeRatePostFlash(:));
sdSpikeRate = std(avgSpikeRatePostFlash(:));
maxSpikeRate = max(avgSpikeRatePostFlash(:));
baselineSpikeRate = mean(avgSpikeRateBaseline(:));

% start with all orientations together
heatAx = axes('Position', [heatLeft btm heatW heatH]); 
hold on;
plotRfMapSmooth(rfmapAllSmoothByOriAll, rfmapAllOriCount, numPixelsPerDegree, mapScale, mapDim, mapXOffset, mapYOffset);
title(sprintf('VES Heatmap (%d-%d ms after flash onset; all orientations)', round(analysisWindowOffset * 1000)), 'Interpreter', 'none');
textParams = {'Units', 'normalized', 'FontSize', 8, 'Color', [1 1 1]};
text(0.02, 0.08, {sprintf('Threshold: %s', thresholdName), ...
        sprintf('Baseline rate: %0.2f Hz', baselineSpikeRate), ...
        sprintf('Mean rate: %0.2f Hz', meanSpikeRate), ...
        sprintf('Std rate: %0.2f Hz', sdSpikeRate), ...
        sprintf('Max rate: %0.2f Hz (%0.1f SDs)', ...
        maxSpikeRate, (maxSpikeRate - meanSpikeRate) / sdSpikeRate)}, ...
        textParams{:});
axis(heatAx, 'square'); % shouldn't do much if i set the dims properly
colorbar;
caxis([0 6]); 

%%
plotFileName = sprintf('%s/%s-%s.png', processedDataDir, unitName, blockName);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
close;
